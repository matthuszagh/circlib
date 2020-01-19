from bisect import insort_left, bisect_left
import concurrent.futures as futures
from subprocess import run
import tempfile
import os
from shutil import rmtree
from typing import List, Tuple
from CSXCAD import CSXCAD as csxcad
import openEMS as openems
from openEMS.physical_constants import C0
from openEMS.ports import Port
import numpy as np


class AutoMesh:
    """
    Automatic mesh generation for OpenEMS CSX structures.

    TODO check and correct if probes in PML.
    TODO ensure at least one mesh line at intersection of metal and probe.
    TODO mesh is subtly wrong (e.g. asymmetrical for cases where it should be
    symmetrical).
    """

    def __init__(
        self,
        csx,
        lmin,
        mres=1 / 20,
        sres=1 / 10,
        smooth=1.5,
        unit=1e-3,
        min_lines=10,
        expand_bounds=[20, 20, 20, 20, 20, 20],
    ):
        """
        @csx is the CSXCAD structure (return value of
        CSXCAD.ContinuousStructure()).
        @lmin is the minimum wavelength associated with the expected
        frequency.
        @mres is the metal resolution, specified as a factor of @lmin.
        @sres is the substrate resolution, specified as a factor of @lmin.
        @smooth is the factor by which adjacent cells are allowed to differ in
        size.
        @unit is the mesh size unit, which defaults to mm.
        @min_lines is the minimum number of mesh lines for a primitive's
        dimensional length, unless the length of that primitive's dimension
        is precisely 0.
        """
        self.csx = csx
        self.lmin = lmin
        self.mres = mres * self.lmin
        self.sres = sres * self.lmin
        self.smooth = smooth
        self.min_lines = min_lines
        # list of 6 elements corresponding to
        # [xmin, xmax, ymin, ymax, zmin, zmax]
        # each element gives the number of cells to add to the mesh
        # at that boundary. The cell size is determined by sres and
        # the actual number of cells added may be more (or possibly)
        # less than what is specified here due to the thirds rule
        # and smoothing. This essentially defines the simulation box.
        # It's anticipated that the user will only define physical
        # structures (e.g. metal layers, substrate, etc.) and will
        # use this to set the simulation box.
        self.expand_bounds = expand_bounds
        # Sort primitives by decreasing priority.
        self.prims = self.csx.GetAllPrimitives()
        # Keep track of mesh regions already applied. This is an array
        # of 3 elements. The 1st element is a list of ranges in the
        # x-dimension that have already been meshed. The 2nd
        # corresponds to the y-dimension and the 3rd corresponds to
        # the z-dimension.
        self.ranges_meshed = [[], [], []]
        # Keep a list of all metal boundaries. No mesh line is allowed to
        # lie on a boundary, and to the extent possible, should obey the
        # thirds rule about that boundary. Zero-dimension metal structures
        # are not added to this list. The 1st item is the x-dimension
        # list, the 2nd is the y-dimension list and the 3rd is the
        # z-dimension list.
        self.metal_bounds = [[], [], []]
        # Mesh lines that cannot be moved. These are mesh lines that
        # lie directly on a zero-dimension primitive.
        self.const_meshes = [[], [], []]
        # Keep track of the smallest valid resolution value. This
        # allows us to later remove all adjacent mesh lines separated
        # by less than this value.
        self.smallest_res = self.mres
        # The generated mesh.
        self.mesh = self.csx.GetGrid()
        self.mesh.SetDeltaUnit(unit)
        # Set the lines first and draw them last since the API doesn't
        # appear to expose a way to remove individual lines.
        self.mesh_lines = [[], [], []]

    def AutoGenMesh(self):
        """
        Start by assuming only two different mesh resolutions: metal and
        substrate/air. This simplifies the 2/3 rule, where the distance to
        the mesh is 1/3 on the metal side and 2/3 on the other side.

        Nonmetal primitives use a mesh resolution of lmin/10 and metal
        primitives use a mesh resolution of lmin/20. There are two
        exceptions to this: (1) if a dimension of a meterial has
        length 0 (e.g. a planar metal sheet) a single mesh line is
        placed exactly on that line, and (2) a nonzero length material
        must have a minimum of 10 mesh lines . The 2nd exception will
        not violate the third's rule or smoothness (that adjacent mesh
        lines not differ in separation by more than a factor of 1.5).
        """
        # add metal mesh
        for prim in self.prims:
            if self._PrimMetalp(prim):
                bounds = self._GetPrimBounds(prim)
                for i in range(3):
                    self._GenMeshInBounds(
                        bounds[i][0], bounds[i][1], self.mres, i, metal=True
                    )

        # add substrate mesh
        for prim in self.prims:
            if not self._PrimMetalp(prim):
                bounds = self._GetPrimBounds(prim)
                for i in range(3):
                    self._GenMeshInBounds(
                        bounds[i][0], bounds[i][1], self.sres, i, metal=False
                    )

        # add simulation box mesh
        for i in range(3):
            self._GenMeshInBounds(
                self.mesh_lines[i][0]
                - (self.sres * self.expand_bounds[2 * i]),
                self.mesh_lines[i][-1]
                + (self.sres * self.expand_bounds[2 * i + 1]),
                self.sres,
                i,
                metal=False,
            )
        # remove unintended, tightly spaced meshes
        for dim in range(3):
            self._RemoveTightMeshLines(dim)

        # enforce thirds rule
        for dim in range(3):
            self._EnforceThirds(dim)

        # smooth mesh
        for dim in range(3):
            self._SmoothMeshLines(dim)

        # set calculated mesh lines
        self._AddAllMeshLines()

    # def ExpandMeshForBoundary(self):
    #     """
    #     Add cells for each boundary where PML is used. The mesh must
    #     already be generated for this to work. @boundary_conds is a
    #     6-element list with the format:

    #     [xmin, xmax, ymin, ymax, zmin, zmax].
    #     """
    #     for boundary, boundary_string in enumerate(self.boundary_conds):
    #         if boundary_string[0:3] == "PML":
    #             num_cells = int(boundary_string[4:])
    #             cell_width = self._GetCellSizeAtBoundary(boundary)
    #             self._AddCellsToBoundary(boundary, num_cells, cell_width)

    # def _GetCellSizeAtBoundary(self, boundary):
    #     """
    #     Compute the cell size at a boundary, where @boundary is an integer
    #     0-5, corresponding to [xmin, xmax, ymin, ymax, zmin, zmax].
    #     """
    #     if self._IsLowerBoundary(boundary):
    #         lower = self.mesh_lines[self._BoundaryDim(boundary)][0]
    #         upper = self.mesh_lines[self._BoundaryDim(boundary)][1]
    #     else:
    #         lower = self.mesh_lines[self._BoundaryDim(boundary)][-2]
    #         upper = self.mesh_lines[self._BoundaryDim(boundary)][-1]
    #     return upper - lower

    # def _AddCellsToBoundary(self, boundary, num_cells, cell_width):
    #     """
    #     Add @num_cells mesh lines with spacing @cell_width to
    #     @boundary. @boundary is an integer from 0-5, corresponding to
    #     [xmin, xmax, ymin, ymax, zmin, zmax].
    #     """
    #     if self._IsLowerBoundary(boundary):
    #         for cell in range(num_cells):
    #             self.mesh_lines[self._BoundaryDim(boundary)].insert(
    #                 0,
    #                 self.mesh_lines[self._BoundaryDim(boundary)][0]
    #                 - cell_width,
    #             )
    #     else:
    #         for cell in range(num_cells):
    #             self.mesh_lines[self._BoundaryDim(boundary)].append(
    #                 self.mesh_lines[self._BoundaryDim(boundary)][-1]
    #                 + cell_width,
    #             )

    # def _IsLowerBoundary(self, boundary):
    #     return boundary == 0 or boundary == 2 or boundary == 4

    # def _BoundaryDim(self, boundary):
    #     return int(boundary / 2)

    # TODO should ensure that inserted mesh lines are not at metal boundaries
    def _EnforceThirds(self, dim):
        for i, pos in enumerate(self.mesh_lines[dim]):
            if (
                pos in self.metal_bounds[dim]
                and pos not in self.const_meshes[dim]
            ):
                # at lower boundary
                if i == 0:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos + (self.mres / 3))
                    self._EnforceThirds(dim)
                # at upper boundary
                elif i == len(self.mesh_lines[dim]) - 1:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos - (self.mres / 3))
                    self._EnforceThirds(dim)
                else:
                    spacing_left = pos - self.mesh_lines[dim][i - 1]
                    spacing_right = self.mesh_lines[dim][i + 1] - pos
                    del self.mesh_lines[dim][i]
                    # metal-metal boundary
                    if spacing_left == spacing_right:
                        new_low = pos - (self.mres / 3)
                        new_high = pos + (self.mres / 3)
                        # no need to add new lines if we already have closer ones.
                        if new_low > pos - spacing_left:
                            insort_left(self.mesh_lines[dim], new_low)
                        if new_high < pos + spacing_right:
                            insort_left(self.mesh_lines[dim], new_high)
                    elif spacing_left < spacing_right:
                        new_low = pos - (self.mres / 3)
                        new_high = pos + (2 * self.mres / 3)
                        if new_low > pos - spacing_left:
                            insort_left(self.mesh_lines[dim], new_low)
                        if new_high < pos + spacing_right:
                            insort_left(self.mesh_lines[dim], new_high)
                    else:
                        new_low = pos - (2 * self.mres / 3)
                        new_high = pos + (self.mres / 3)
                        if new_low > pos - spacing_left:
                            insort_left(self.mesh_lines[dim], new_low)
                        if new_high < pos + spacing_right:
                            insort_left(self.mesh_lines[dim], new_high)
                    self._EnforceThirds(dim)

    def _RemoveTightMeshLines(self, dim):
        """
        Remove adjacent mesh lines for dimension @dim with spacing less
        than the smallest valid resolution.
        """
        last_pos = self.mesh_lines[dim][0]
        for i, pos in enumerate(self.mesh_lines[dim]):
            if i == 0:
                continue
            # we can freely delete duplicates
            if pos == last_pos:
                del self.mesh_lines[dim][i]
                self._RemoveTightMeshLines(dim)
            # we have to check whether these are zero-dimension
            # structures before deleting them.
            elif pos - last_pos < self.smallest_res and (
                pos not in self.const_meshes[dim]
                or last_pos not in self.const_meshes[dim]
            ):
                if last_pos not in self.const_meshes[dim]:
                    del self.mesh_lines[dim][i - 1]
                else:
                    del self.mesh_lines[dim][i]
                self._RemoveTightMeshLines(dim)
            else:
                last_pos = pos

    def _AddAllMeshLines(self):
        for dim in range(3):
            for line in self.mesh_lines[dim]:
                self.mesh.AddLine(dim, line)

    def _GetMesh(self):
        return self.mesh

    def _PrimMetalp(self, prim):
        return prim.GetProperty().GetTypeString() == "Metal"

    def _GetPrimBounds(self, prim):
        orig_bounds = prim.GetBoundBox()
        bounds = [[None, None], [None, None], [None, None]]
        for i in range(3):
            upper = max(orig_bounds[0][i], orig_bounds[1][i])
            lower = min(orig_bounds[0][i], orig_bounds[1][i])
            bounds[i] = [lower, upper]
        return bounds

    def _MeshResInBounds(self, lower, upper, dim):
        lower_idx = bisect_left(self.mesh_lines[dim], lower)
        upper_idx = min(
            bisect_left(self.mesh_lines[dim], upper) + 1,
            len(self.mesh_lines[dim]),
        )
        spacing = []
        last_pos = self.mesh_lines[dim][lower_idx]
        for idx in range(lower_idx + 1, upper_idx):
            spacing.append(self.mesh_lines[dim][idx] - last_pos)
            last_pos = self.mesh_lines[dim][idx]
        return sum(spacing) / len(spacing)

    def _SplitBounds(self, lower, upper, dim):
        """
        Split bounds delimited by [lower, upper] into regions where mesh
        already exists and regions where it doesn't yet exist.

        Returns a list of 2 items. The 1st item is a list of bounds
        where the new mesh boundaries are outside the existing
        mesh. The 2nd item is a list of bounds where the new mesh
        overlaps the existing mesh.
        """
        outin_ranges = [[], []]
        if len(self.ranges_meshed[dim]) == 0:
            outin_ranges[0].append([lower, upper])
            return outin_ranges

        for [lower_mesh, upper_mesh] in self.ranges_meshed[dim]:
            if upper <= lower_mesh:
                outin_ranges[0].append([lower, upper])
                # since meshed ranges are sorted, we can ignore the rest.
                return outin_ranges
            elif lower >= upper_mesh:
                continue
            elif lower < lower_mesh:
                outin_ranges[0].append([lower, lower_mesh])
                if upper > upper_mesh:
                    outin_ranges[1].append([lower_mesh, upper_mesh])
                    lower = upper_mesh
                    continue
                else:
                    outin_ranges[1].append([lower_mesh, upper])
                    return outin_ranges
            else:
                outin_ranges[1].append([lower, min(upper, upper_mesh)])
                if upper > upper_mesh:
                    lower = upper_mesh
                    continue
                else:
                    return outin_ranges
        if lower < upper:
            outin_ranges[0].append([lower, upper])

        return outin_ranges

    def _ClearMeshInBounds(self, lower, upper, dim):
        for elt in self.mesh_lines[dim]:
            if elt >= lower and elt <= upper:
                self.mesh_lines[dim].remove(elt)

    def _RangeUnion(self, ranges, start_idx=0):
        ranges = sorted(ranges)
        if len(ranges[start_idx:]) <= 1:
            return ranges

        if ranges[start_idx][1] >= ranges[start_idx + 1][0]:
            ranges.append([ranges[start_idx][0], ranges[start_idx + 1][1]])
            del ranges[start_idx : start_idx + 2]
            return self._RangeUnion(ranges[start_idx:])
        else:
            return self._RangeUnion(ranges[start_idx + 1 :])

    def _ConsolidateMeshedRanges(self, dim):
        """
        Order meshed ranges and consolidate contiguous ranges.
        """
        self.ranges_meshed[dim] = self._RangeUnion(self.ranges_meshed[dim])

    def _UpdateRanges(self, lower, upper, dim):
        """
        @dim is the dimension: 0, 1, 2 for x, y, or z.
        """
        self.ranges_meshed[dim].append([lower, upper])
        self._ConsolidateMeshedRanges(dim)

    def _SmoothMeshLines(self, dim):
        for i, pos in enumerate(self.mesh_lines[dim]):
            if i == 0 or i == len(self.mesh_lines[dim]) - 1:
                continue
            left_spacing = pos - self.mesh_lines[dim][i - 1]
            right_spacing = self.mesh_lines[dim][i + 1] - pos
            if (
                left_spacing > (self.smooth * right_spacing)
                and left_spacing - (self.smooth * right_spacing)
                > right_spacing / 10
            ):
                if left_spacing / right_spacing <= 2:
                    adj = (left_spacing - (self.smooth * right_spacing)) / (
                        self.smooth + 1
                    )
                    # TODO need to ensure new mesh line doesn't fall
                    # on metal boundary or violate thirds.
                    if pos not in self.const_meshes[dim]:
                        del self.mesh_lines[dim][i]
                        insort_left(self.mesh_lines[dim], pos - adj)
                    else:
                        insort_left(
                            self.mesh_lines[dim], pos - (left_spacing / 2)
                        )
                else:
                    insort_left(self.mesh_lines[dim], pos - (left_spacing / 2))
                self._SmoothMeshLines(dim)
            elif (
                right_spacing > self.smooth * left_spacing
                and right_spacing - (self.smooth * left_spacing)
                > left_spacing / 10
            ):
                if right_spacing / left_spacing <= 2:
                    adj = (right_spacing - (self.smooth * left_spacing)) / (
                        self.smooth + 1
                    )
                    # TODO need to ensure new mesh line doesn't fall
                    # on metal boundary or violate thirds.
                    if pos not in self.const_meshes[dim]:
                        del self.mesh_lines[dim][i]
                        insort_left(self.mesh_lines[dim], pos + adj)
                    else:
                        insort_left(
                            self.mesh_lines[dim], pos + (right_spacing / 2)
                        )
                else:
                    insort_left(
                        self.mesh_lines[dim], pos + (right_spacing / 2)
                    )
                self._SmoothMeshLines(dim)

    def _NearestDivisibleRes(self, lower, upper, res):
        """
        Return the nearest resolution to @res that evenly subdivides the
        interval [@lower, @upper].

        This is important because it helps prevent adjacent lines from
        bunching up and unnecessarily increasing the simulation time.
        """
        num_divisions = np.round((upper - lower) / res)
        num_divisions = max(num_divisions, 1)
        return (upper - lower) / num_divisions

    def _GenMeshInBounds(self, lower, upper, res, dim, metal=False):
        if lower == upper:
            insort_left(self.mesh_lines[dim], lower)
            insort_left(self.const_meshes[dim], lower)
        else:
            [outer_bounds, inner_bounds] = self._SplitBounds(lower, upper, dim)
            for obound in outer_bounds:
                if upper - lower < self.min_lines * res:
                    res = (upper - lower) / self.min_lines
                else:
                    res = self._NearestDivisibleRes(obound[0], obound[1], res)
                self.smallest_res = min(self.smallest_res, res)
                j = obound[0]
                while j <= obound[1]:
                    insort_left(self.mesh_lines[dim], j)
                    j += res
                self._UpdateRanges(obound[0], obound[1], dim)
                if metal:
                    insort_left(self.metal_bounds[dim], obound[0])
                    insort_left(self.metal_bounds[dim], obound[1])
            for ibound in inner_bounds:
                if upper - lower < self.min_lines * res:
                    res = (upper - lower) / self.min_lines
                else:
                    res = self._NearestDivisibleRes(ibound[0], ibound[1], res)
                self.smallest_res = min(self.smallest_res, res)
                # only redo the mesh if the desired one is finer than
                # the existing one
                cur_mesh_res = self._MeshResInBounds(ibound[0], ibound[1], dim)
                if cur_mesh_res > res and abs(cur_mesh_res - res) > res / 10:
                    self._ClearMeshInBounds(ibound[0], ibound[1], dim)
                    j = ibound[0]
                    while j <= ibound[1]:
                        insort_left(self.mesh_lines[dim], j)
                        j += res
                    self._UpdateRanges(ibound[0], ibound[1], dim)
                if metal:
                    insort_left(self.metal_bounds[dim], ibound[0])
                    insort_left(self.metal_bounds[dim], ibound[1])


class PCB:
    """
    A PCB structure and material properties for use in an openems simulation.
    """

    def __init__(
        self,
        layers: int,
        sub_epsr: List[Tuple[float, float]],
        sub_rho: float,
        layer_sep: List[float],
        layer_thickness: List[float],
    ):
        """
        Args:
            layers: number of conductive layers
            sub_epsr: substrate dielectric constant. dictionary of frequency
                      (Hz) and associated dielectric.
            sub_rho: volume resistivity (ohm*mm)
            layer_sep: separations (in mm) between adjacent copper layers. A
                       list where the first value is the separation between
                       the top layer and second layer, etc. This is
                       equivalently the substrate thickness.
            layer_thickness: thickness of each conductive layer (in mm). Again
                             proceeds from top to bottom layer.
        """
        self.layers = layers
        self.sub_epsr = sorted(sub_epsr, key=lambda epsr: epsr[0])
        self.sub_rho = sub_rho
        self.sub_kappa = 1 / sub_rho
        self.layer_sep = layer_sep
        self.layer_thickness = layer_thickness

    def epsr_at_freq(self, freq: float):
        """
        Approximate the dielectric at a given frequency given the provided
        epsr values.

        Args:
            freq: frequency of interest (Hz)

        Returns:
            dielectric constant
        """
        if freq <= self.sub_epsr[0][0]:
            return self.sub_epsr[0][1]
        elif freq >= self.sub_epsr[-1][0]:
            return self.sub_epsr[-1][1]

        # perform linear interpolation
        tup_low = bisect_left([x[0] for x in self.sub_epsr], freq)
        if self.sub_epsr[tup_low][0] == freq:
            return self.sub_epsr[tup_low][1]

        tup_high = tup_low
        tup_low -= 1
        xlow = self.sub_epsr[tup_low][0]
        xhigh = self.sub_epsr[tup_high][0]
        ylow = self.sub_epsr[tup_low][1]
        yhigh = self.sub_epsr[tup_high][1]
        slope = (yhigh - ylow) / (xhigh - xlow)
        return ylow + (slope * (freq - xlow))


common_pcbs = {
    "oshpark4": PCB(
        layers=4,
        sub_epsr=[
            (100e6, 3.72),
            (1e9, 3.69),
            (2e9, 3.68),
            (5e9, 3.64),
            (10e9, 3.65),
        ],
        sub_rho=4.4e14,
        layer_sep=[0.1702, 1.1938, 0.1702],
        layer_thickness=[0.0356, 0.0178, 0.0178, 0.0356],
    )
}


def wheeler_z0(w: float, t: float, er: float, h: float) -> float:
    """
    Calculate the microstrip characteristic impedance for a given width
    using Wheeler's equation. Wheeler's equation can be found at:

    https://en.wikipedia.org/wiki/Microstrip#Characteristic_impedance

    Args:
        w: microstrip trace width (mm)
        t: trace thickness (mm)
        er: substrate relative permittivity
        h: substrate height (thickness) (mm)

    Returns:
        characteristic impedance
    """
    z0 = 376.730313668
    w *= 1e-3
    t *= 1e-3
    h *= 1e-3
    weff = w + (
        (t * ((1 + (1 / er)) / (2 * np.pi)))
        * np.log(
            (4 * np.e)
            / (
                np.sqrt(
                    ((t / h) ** 2)
                    + (((1 / np.pi) * ((1 / ((w / t) + (11 / 10))))) ** 2)
                )
            )
        )
    )
    tmp1 = 4 * h / weff
    tmp2 = (14 + (8 / er)) / 11
    zm = (z0 / (2 * np.pi * np.sqrt(2 * (1 + er)))) * np.log(
        1
        + (
            tmp1
            * (
                (tmp2 * tmp1)
                + (
                    np.sqrt(
                        (
                            (tmp2 * tmp1) ** 2
                            + ((np.pi ** 2) * ((1 + (1 / er)) / 2))
                        )
                    )
                )
            )
        )
    )
    return zm


def wheeler_z0_width(
    z0: float,
    t: float,
    er: float,
    h: float,
    tol: float = 0.01,
    guess: float = 0.3,
) -> float:
    """
    Calculate the microstrip width for a given characteristic
    impedance using Wheeler's formula.

    Args:
        z0: characteristic impedance (ohm)
        t: trace thickness (mm)
        er: substrate relative permittivity
        h: substrate height (thickness) (mm)
        tol: acceptable impedance tolerance (ohm)
        guess: an initial guess for the width (mm). This can improve
               convergence time when the approximate width is known.

    Returns:
        trace width (mm)
    """
    width = guess
    zm = wheeler_z0(w=width, t=t, er=er, h=h)
    wlow = width / 10
    zlow = wheeler_z0(w=wlow, t=t, er=er, h=h)
    # inverse relation between width and z0
    while zlow < z0:
        wlow /= 10
        zlow = wheeler_z0(w=wlow, t=t, er=er, h=h)

    whigh = width * 10
    zhigh = wheeler_z0(w=whigh, t=t, er=er, h=h)
    while zhigh > z0:
        whigh *= 10
        zhigh = wheeler_z0(w=whigh, t=t, er=er, h=h)

    while np.absolute(zm - z0) > tol:
        if zm > z0:
            m = (zhigh - zm) / (whigh - width)
            wlow = width
            zlow = zm
        else:
            m = (zm - zlow) / (width - wlow)
            whigh = width
            zhigh = zm

        # use linear interpolation to update guess
        width = width + ((z0 - zm) / m)
        zm = wheeler_z0(w=width, t=t, er=er, h=h)

    return width


class Microstrip:
    """
    Microstrip transmission line.
    """

    def __init__(
        self,
        pcb: PCB,
        f0: float,
        fc: float,
        z0_ref: float = 50,
        microstrip_width: float = None,
        microstrip_len: float = 100,
        substrate_width: float = 20,
        substrate_len: float = 120,
        fcsx: str = None,
        fvtr_dir: str = None,
        efield: bool = False,
    ):
        """
        Args:
            pcb: PCB object
            f0: center frequency (Hz)
            fc: corner frequency, given as difference from f0 (Hz)
            z0_ref: desired impedance (ohm)
            microstrip_width: microstrip trace width (mm). If ommitted, an
                              analytical best guess for z0_ref and f0 will be
                              used (see `wheeler_z0_width`).
            microstrip_len: microstrip trace length (mm). This value isn't
                            critical but does need to be large enough for the
                            openems simulation to work. If unsure, stick with
                            the default.
            substrate_width: substrate width (mm)
            substrate_len: substrate length (mm). Must be at least as large as
                           the microstrip length. If unsure, use the default.
            fcsx: CSX file name to write CSXCAD structure to. If ommitted, do
                  not save the file. Omitting this does not prevent you from
                  viewing the CSX structure, just from saving it for later.
                  Relative to current directory.
            fvtr_dir: directory for VTR E-field dump files. Files will be
                      prefixed with 'Et_'. If ommitted, E-field can still be
                      viewed unless `efield` is set to false in which case this
                      parameter has no effect.
            efield: dump E-field time values for viewing.
        """
        self.pcb = pcb
        self.f0 = f0
        self.fc = fc
        self.z0_ref = z0_ref
        if microstrip_width is None:
            self.microstrip_width = wheeler_z0_width(
                z0=z0_ref,
                t=self.pcb.layer_thickness[0],
                er=self.pcb.epsr_at_freq(f0),
                h=self.pcb.layer_sep[0],
            )
        else:
            self.microstrip_width = microstrip_width
        self.microstrip_len = microstrip_len
        self.substrate_width = substrate_width
        self.substrate_len = substrate_len
        if fcsx is None:
            fcsx = tempfile.TemporaryFile().name
        else:
            self.fcsx = os.path.abspath(fcsx)

        self.efield = efield
        if fvtr_dir is None:
            if efield is not None:
                fvtr_dir = tempfile.TemporaryDirectory().name
        else:
            self.fvtr_dir = os.path.abspath(fvtr_dir)
            if os.path.exists(self.fvtr_dir):
                rmtree(self.fvtr_dir)
            os.mkdir(self.fvtr_dir)

        # track whether simulation already performed so that we only
        # need to simulate once
        self.sim_done = False
        self.csx_done = False
        self.freq_bins = None
        self.calc_ports = None
        self.fdtd = None

    def sim(
        self,
        num_probes: int = 3,
        num_freq_bins: int = 500,
        zero_trace_height: bool = False,
    ) -> None:
        """
        Run an openems simulation for the current microstrip structure and
        write the results to a file.

        Args:
            num_probes: number of impedance probes placed along trace. This
                        also determines outer dimension of return value.
            num_freq_bins: number of frequency bins. More frequency bins means
                       longer simulation time but possibly greater accuracy.
            zero_trace_height: approximate trace height as 0. This can greatly
                               reduce simulation time but may reduce accuracy.
        """
        self.gen_csx(
            num_probes=num_probes, zero_trace_height=zero_trace_height
        )
        tmpdir = tempfile.TemporaryDirectory()
        self.fdtd.Run(tmpdir.name, cleanup=True)

        self.freq_bins = np.linspace(
            self.f0 - self.fc, self.f0 + self.fc, num_freq_bins
        )
        self.calc_ports = [None] * num_probes
        for i, _ in enumerate(self.calc_ports):
            self.calc_ports[i] = Port(
                None,
                None,
                None,
                None,
                None,
                U_filenames=["ut_" + str(i)],
                I_filenames=["it_" + str(i)],
            )
            self.calc_ports[i].CalcPort(
                tmpdir.name, self.freq_bins, ref_impedance=self.z0_ref
            )

        self.sim_done = True

    def gen_csx(
        self, num_probes: int = 3, zero_trace_height: bool = False,
    ) -> None:
        """
        Generate CSX structure.

        Args:
            num_probes: number of impedance probes placed along trace. This
                        also determines outer dimension of return value.
            zero_trace_height: approximate trace height as 0. This can greatly
                               reduce simulation time but may reduce accuracy.
        """
        # dimensions and parameters
        unit = 1e-3
        lmin = C0 / ((self.f0 + self.fc) * unit)
        if zero_trace_height:
            trace_height = 0
        else:
            trace_height = self.pcb.layer_thickness[0]

        # structures
        csx = csxcad.ContinuousStructure()
        microstrip = csx.AddMetal("PEC")
        microstrip.AddBox(
            priority=10,
            start=[-self.microstrip_len / 2, -self.microstrip_width / 2, 0],
            stop=[
                self.microstrip_len / 2,
                self.microstrip_width / 2,
                trace_height,
            ],
        )
        substrate = csx.AddMaterial(
            "substrate",
            epsilon=self.pcb.epsr_at_freq(self.f0),
            kappa=self.pcb.sub_kappa,
        )
        substrate.AddBox(
            priority=0,
            start=[
                -self.substrate_len / 2,
                -self.substrate_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[self.microstrip_len / 2, self.substrate_width / 2, 0],
        )

        # simulation
        self.fdtd = openems.openEMS(EndCriteria=1e-5)
        self.fdtd.SetCSX(csx)
        self.fdtd.SetGaussExcite(self.f0, self.fc)
        self.fdtd.SetBoundaryCond(
            ["PML_8", "PML_8", "MUR", "MUR", "PEC", "MUR"]
        )
        self.fdtd.AddLumpedPort(
            port_nr=0,
            R=self.z0_ref,
            start=[
                -self.microstrip_len / 2 + (self.microstrip_len / 200),
                -self.microstrip_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[
                -self.microstrip_len / 2 + (self.microstrip_len / 200),
                self.microstrip_width / 2,
                0,
            ],
            p_dir="z",
            excite=1,
            priority=999,
        )
        vprobe = [None] * num_probes
        iprobe = [None] * num_probes

        # TODO num_probes here and iprobe doesn't quite work because
        # probes could be placed in PML
        for i, _ in enumerate(vprobe):
            vprobe[i] = csx.AddProbe("ut_" + str(i), 0)
            vprobe[i].AddBox(
                start=[
                    -self.microstrip_len / 2
                    + ((i + 1) * self.microstrip_len / (num_probes + 1)),
                    0,
                    -self.pcb.layer_sep[0],
                ],
                stop=[
                    -self.microstrip_len / 2
                    + ((i + 1) * self.microstrip_len / (num_probes + 1)),
                    0,
                    0,
                ],
            )

        for i, _ in enumerate(iprobe):
            iprobe[i] = csx.AddProbe("it_" + str(i), 1, norm_dir=0)
            iprobe[i].AddBox(
                start=[
                    -self.microstrip_len / 2
                    + ((i + 1) * self.microstrip_len / (num_probes + 1)),
                    -self.microstrip_width / 2,
                    0,
                ],
                stop=[
                    -self.microstrip_len / 2
                    + ((i + 1) * self.microstrip_len / (num_probes + 1)),
                    self.microstrip_width / 2,
                    trace_height,
                ],
            )

        auto_mesh = AutoMesh(
            csx,
            lmin,
            mres=1 / 20,
            sres=1 / 10,
            smooth=1.4,
            unit=unit,
            min_lines=5,
            expand_bounds=[10, 0, 10, 10, 0, 10],
        )
        auto_mesh.AutoGenMesh()

        # E-field recording
        if self.efield:
            Et = csx.AddDump(os.path.join(self.fvtr_dir, "Et_"), file_type=0)
            start = [
                -self.substrate_len / 2,
                -self.substrate_width / 2,
                -self.pcb.layer_sep[0] / 2,
            ]
            stop = [
                self.microstrip_len / 2,
                self.substrate_width / 2,
                -self.pcb.layer_sep[0] / 2,
            ]
            Et.AddBox(start, stop)

        # write CSX file
        csx.Write2XML(self.fcsx)
        self.csx_done = True

    def port_z0(self) -> List[List[float]]:
        """
        Get impedance as a function of frequency for each port.
        """
        if not self.sim_done:
            raise RuntimeWarning(
                "Must run simulation before retreiving impedance values. "
                "Doing nothing."
            )
        else:
            z0 = [
                np.absolute(self.calc_ports[i].uf_tot)
                / np.absolute(self.calc_ports[i].if_tot)
                for i, _ in enumerate(self.calc_ports)
            ]

            return z0

    def avg_port_z0(self) -> List[float]:
        """
        Get impedance as a function of frequency averaged across all ports.
        """
        z0 = self.port_z0()
        return np.average(z0, axis=0)

    def view_csx(self):
        """
        View CSX structure. This blocks the calling thread and AppCSXCAD
        must be closed before proceeding.
        """
        if not self.csx_done:
            raise RuntimeWarning(
                "Must generate CSX before viewing it. Doing nothing."
            )
        else:
            run(["AppCSXCAD", self.fcsx])

    def view_efield(self):
        """
        View E-field time dump. This blocks the calling thread and
        Paraview must be closed before proceeding.
        """
        if not self.sim_done:
            raise RuntimeWarning(
                "Must run simulation before viewing E-field dump. "
                "Doing nothing."
            )
        else:
            run(
                [
                    "paraview",
                    "--data={}".format(
                        os.path.join(self.fvtr_dir, "Et__..vtr")
                    ),
                ]
            )


def microstrip_width_sweep(
    f0: float,
    fc: float,
    thick: float,
    epsr: float,
    resistivity: float,
    z0_ref: float,
    width_bounds: [float, float],
    num_points: int,
    fout: str,
) -> None:
    """
    Calculate microstrip characteristic impedance as a function of trace width.

    Args:
        f0: center frequency (Hz)
        fc: 20dB cutoff frequency, as the deviation from the center frequency
        thick: substrate thickness (mm)
        epsr: substrate dielectric constant
        z0_ref: desired impedance
        resistivity: substrate resistivity
        width_bounds: [lower_width, upper_width] sweep over this width range
        num_points: number of width values to calculate. This value has a
                    significant effect on simulation time since each point
                    is calculated in its own thread. Therefore, to minimize
                    computation time, its recommended to choose some multiple
                    of the number of cores available on the simulation machine.
                    The total simulation time will be roughly equal to the time
                    it takes to compute 1 point times the ratio of the number
                    of points to number of machine cores.
        fout: filename where results should be written. relative to current
              dir.
    """
    widths = np.linspace(width_bounds[0], width_bounds[1], num_points)
    subprocs = [None] * len(widths)
    executer = futures.ProcessPoolExecutor()
    for i, width in enumerate(widths):
        subprocs[i] = executer.submit(
            microstrip_impedance,
            f0,
            fc,
            thick,
            epsr,
            resistivity,
            width,
            z0_ref,
        )

    with open(fout, "w") as f:
        f.write(
            "{:10} {:10} {:10} {:10}\n".format("width", "zlow", "z0", "zhigh")
        )
        for i, _ in enumerate(subprocs):
            [zlow, z0, zhigh] = subprocs[i].result()
            f.write(
                "{:10f} {:10f} {:10f} {:10f}\n".format(
                    widths[i], zlow, z0, zhigh
                )
            )


def microstrip_impedance(
    f0: float,
    fc: float,
    thick: float,
    epsr: float,
    resistivity: float,
    width: float,
    z0_ref: float,
) -> [float, float, float]:
    """
    Calculate the characeristic impedance for a microstrip line for a
    frequency range and given width.

    Returns:
        [z0_flow, z0_fcenter, z0_fhigh] where flow=f0-fc, fcenter=f0,
        and fhigh=f0+fc.
    """
    # dimensions and parameters
    unit = 1e-3
    microstrip_len = 100
    substrate_width = 20
    substrate_len = 1.2 * microstrip_len
    lmin = C0 / ((f0 + fc) * unit)
    substrate_kappa = 1 / resistivity

    # structures
    csx = csxcad.ContinuousStructure()
    microstrip = csx.AddMetal("PEC")
    microstrip.AddBox(
        priority=10,
        start=[-microstrip_len / 2, -width / 2, 0],
        stop=[microstrip_len / 2, width / 2, 0],
    )
    substrate = csx.AddMaterial(
        "substrate", epsilon=epsr, kappa=substrate_kappa
    )
    substrate.AddBox(
        priority=0,
        start=[-substrate_len / 2, -substrate_width / 2, -thick],
        stop=[substrate_len / 2, substrate_width / 2, 0],
    )

    # simulation
    fdtd = openems.openEMS(EndCriteria=1e-5)
    fdtd.SetCSX(csx)
    fdtd.SetGaussExcite(f0, fc)
    fdtd.SetBoundaryCond(["PML_8", "PML_8", "MUR", "MUR", "PEC", "MUR"])
    fdtd.AddLumpedPort(
        port_nr=0,
        R=z0_ref,
        start=[-microstrip_len / 2 + (microstrip_len / 200), 0, -thick],
        stop=[-microstrip_len / 2 + (microstrip_len / 200), 0, 0],
        p_dir="z",
        excite=1,
        priority=999,
    )
    vprobe = [None, None, None]
    iprobe = [None, None, None]

    for i, probe in enumerate(vprobe):
        vprobe[i] = csx.AddProbe("vprobe_" + str(width) + "_" + str(i), 0)
        vprobe[i].AddBox(
            start=[
                -microstrip_len / 2 + ((i + 1) * microstrip_len / 4),
                0,
                -thick,
            ],
            stop=[-microstrip_len / 2 + ((i + 1) * microstrip_len / 4), 0, 0,],
        )

    for i, probe in enumerate(iprobe):
        iprobe[i] = csx.AddProbe(
            "iprobe_" + str(width) + "_" + str(i), 1, norm_dir=0
        )
        iprobe[i].AddBox(
            start=[
                -microstrip_len / 2
                - (microstrip_len / 20)
                + ((i + 1) * microstrip_len / 4),
                -width / 2,
                0,
            ],
            stop=[
                -microstrip_len / 2
                - (microstrip_len / 20)
                + ((i + 1) * microstrip_len / 4),
                width / 2,
                0,
            ],
        )

    auto_mesh = AutoMesh(
        csx,
        lmin,
        mres=1 / 20,
        sres=1 / 10,
        smooth=1.4,
        unit=unit,
        min_lines=5,
        expand_bounds=[10, 0, 10, 10, 0, 10],
    )
    auto_mesh.AutoGenMesh()

    tmpdir = tempfile.TemporaryDirectory()
    fdtd.Run(tmpdir.name, cleanup=True)

    num_freq = 2000
    freq = np.linspace(f0 - fc, f0 + fc, num_freq)
    z0_ports = [None, None, None]
    for n, _ in enumerate(z0_ports):
        z0_ports[n] = Port(
            None,
            None,
            None,
            None,
            None,
            U_filenames=["vprobe_" + str(width) + "_" + str(n)],
            I_filenames=["iprobe_" + str(width) + "_" + str(n)],
        )
        z0_ports[n].CalcPort(tmpdir.name, freq, ref_impedance=z0_ref)

    zlow = np.average(
        [
            np.absolute(z0_ports[i].uf_tot[0])
            / np.absolute(z0_ports[i].if_tot[0])
            for i in range(3)
        ]
    )
    z0 = np.average(
        [
            np.absolute(z0_ports[i].uf_tot[int(num_freq / 2 - 1)])
            / np.absolute(z0_ports[i].if_tot[int(num_freq / 2 - 1)])
            for i in range(3)
        ]
    )
    zhigh = np.average(
        [
            np.absolute(z0_ports[i].uf_tot[-1])
            / np.absolute(z0_ports[i].if_tot[-1])
            for i in range(3)
        ]
    )
    return [zlow, z0, zhigh]
