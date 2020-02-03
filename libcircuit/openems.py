from bisect import insort_left, bisect_left
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial
from subprocess import run
import tempfile
import os
from shutil import rmtree
from typing import List, Tuple
from CSXCAD import CSXCAD as csxcad
import openEMS as openems
from openEMS.physical_constants import C0
from openEMS.ports import UI_data
import numpy as np


class Mesh:
    """
    Automatic mesh generation for OpenEMS CSX structures.  Probes
    should always be defined after mesh generation.  Additionally,
    voltage probes should be placed on a mesh line and current probes
    should be placed midway between two adjacent mesh lines.

    TODO mesh is subtly wrong (e.g. asymmetrical for cases where it
    should be symmetrical).
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
        :param csx: the CSXCAD structure (return value of
            CSXCAD.ContinuousStructure()).
        :param lmin: the minimum wavelength associated with the
            expected frequency.
        :param mres: the metal resolution, specified as a factor of
            lmin.
        :param sres: the substrate resolution, specified as a factor
            of lmin.
        :param smooth: the factor by which adjacent cells are allowed
            to differ in size.
        :param unit: the mesh size unit, which defaults to mm.
        :param min_lines: the minimum number of mesh lines for a
            primitive's dimensional length, unless the length of that
            primitive's dimension is precisely 0.
        """
        self.csx = csx
        self.lmin = lmin
        self.mres = mres * self.lmin
        self.sres = sres * self.lmin
        self.smooth = smooth
        # mesh lines are added at both boundaries, which gives us an
        # extra mesh line
        self.min_lines = min_lines - 1
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

    def generate_mesh(self, enforce_thirds=True, smooth=True):
        """
        Start by assuming only two different mesh resolutions: metal
        and substrate/air.  This simplifies the 2/3 rule, where the
        distance to the mesh is 1/3 on the metal side and 2/3 on the
        other side.

        Nonmetal primitives use a mesh resolution of lmin/10 and metal
        primitives use a mesh resolution of lmin/20.  There are two
        exceptions to this: (1) if a dimension of a meterial has
        length 0 (e.g. a planar metal sheet) a single mesh line is
        placed exactly on that line, and (2) a nonzero length material
        must have a minimum of 10 mesh lines .  The 2nd exception will
        not violate the third's rule or smoothness (that adjacent mesh
        lines not differ in separation by more than a factor of 1.5).

        :param enforce_thirds: Enforce thirds rule for metal
            boundaries.  This should always be enabled unless you want
            to debug the mesh generation.
        :param smooth: Smooth mesh lines so that adjacent separations
            do not differ by more than the smoothness factor.  This
            should always be enabled unless you want to debug the mesh
            generation.
        """
        # add metal mesh
        for prim in self.prims:
            if self._type_str(prim) == "Metal":
                bounds = self._get_prim_bounds(prim)
                for i in range(3):
                    self._gen_mesh_in_bounds(
                        bounds[i][0], bounds[i][1], self.mres, i, metal=True
                    )

        # add substrate mesh
        for prim in self.prims:
            if self._type_str(prim) == "Material":
                bounds = self._get_prim_bounds(prim)
                for i in range(3):
                    self._gen_mesh_in_bounds(
                        bounds[i][0], bounds[i][1], self.sres, i, metal=False
                    )

        # add simulation box mesh
        for i in range(3):
            self._gen_mesh_in_bounds(
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
            self._remove_tight_mesh_lines(dim)

        # enforce thirds rule
        if enforce_thirds:
            for dim in range(3):
                self._enforce_thirds(dim)

        # smooth mesh
        if smooth:
            for dim in range(3):
                self._smooth_mesh_lines(dim)

        # set calculated mesh lines
        self._add_lines_to_mesh()

    def nearest_mesh_line(self, dim: int, pos: float) -> (int, float):
        """
        Find the nearest mesh line to a desired position for a given
        dimension.

        :param dim: 0, 1, or 2 for x, y, z.
        :param pos: desired position.

        :returns: (index, position) where index is the array index and
                  position is the actual dimension value.
        """
        lines = self.mesh_lines[dim]
        bisect_pos = bisect_left(self.mesh_lines[dim], pos)
        if bisect_pos == 0:
            return (0, lines[0])
        elif bisect_pos == len(lines):
            return (bisect_pos - 1, lines[bisect_pos - 1])
        else:
            lower = lines[bisect_pos - 1]
            upper = lines[bisect_pos]
            if pos - lower < upper - pos:
                return (bisect_pos - 1, lower)
            else:
                return (bisect_pos, upper)

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
    def _enforce_thirds(self, dim):
        """
        Replace mesh lines at metal boundaries with a mesh line
        1/3*res inside the metal boundary and 2/3*res outside.

        :param dim: Dimension for which thirds rule should be
            enforced.
        """
        for i, pos in enumerate(self.mesh_lines[dim]):
            if (
                pos in self.metal_bounds[dim]
                and pos not in self.const_meshes[dim]
            ):
                # at lower boundary
                if i == 0:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos + (self.mres / 3))
                    self._enforce_thirds(dim)
                # at upper boundary
                elif i == len(self.mesh_lines[dim]) - 1:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos - (self.mres / 3))
                    self._enforce_thirds(dim)
                else:
                    spacing_left = pos - self.mesh_lines[dim][i - 1]
                    spacing_right = self.mesh_lines[dim][i + 1] - pos
                    del self.mesh_lines[dim][i]
                    # metal-metal boundary
                    if (
                        abs(spacing_left - spacing_right)
                        < self.smallest_res / 10
                    ):
                        new_low = pos - (spacing_left / 2)
                        new_high = pos + (spacing_left / 2)
                    # don't need to add tolerance for float comparison
                    # since metal-metal boundary check already did
                    # that
                    elif spacing_left < spacing_right:
                        new_low = pos - (spacing_left / 3)
                        new_high = pos + (2 * spacing_left / 3)
                    else:
                        new_low = pos - (2 * spacing_right / 3)
                        new_high = pos + (spacing_right / 3)

                    insort_left(self.mesh_lines[dim], new_low)
                    insort_left(self.mesh_lines[dim], new_high)
                    self._enforce_thirds(dim)

    def _remove_tight_mesh_lines(self, dim):
        """
        Remove adjacent mesh lines for dimension @dim with spacing
        less than the smallest valid resolution.

        :param dim: Dimension in which to remove tightly spaced
            meshes.
        """
        last_pos = self.mesh_lines[dim][0]
        for i, pos in enumerate(self.mesh_lines[dim]):
            if i == 0:
                continue
            # we can freely delete duplicates
            if pos == last_pos:
                del self.mesh_lines[dim][i]
                self._remove_tight_mesh_lines(dim)
            # we have to check whether these are zero-dimension
            # structures before deleting them.
            elif (
                pos - last_pos < self.smallest_res
                and abs(pos - last_pos - self.smallest_res)
                > self.smallest_res / 10
                and (
                    pos not in self.const_meshes[dim]
                    or last_pos not in self.const_meshes[dim]
                )
            ):
                if last_pos not in self.const_meshes[dim]:
                    del self.mesh_lines[dim][i - 1]
                else:
                    del self.mesh_lines[dim][i]
                self._remove_tight_mesh_lines(dim)
            else:
                last_pos = pos

    def _add_lines_to_mesh(self):
        """
        Generates the actual CSX mesh structure from lines.
        """
        for dim in range(3):
            for line in self.mesh_lines[dim]:
                self.mesh.AddLine(dim, line)

    def _get_mesh(self):
        return self.mesh

    def _type_str(self, prim):
        return prim.GetProperty().GetTypeString()

    def _get_prim_bounds(self, prim):
        orig_bounds = prim.GetBoundBox()
        bounds = [[None, None], [None, None], [None, None]]
        for i in range(3):
            upper = max(orig_bounds[0][i], orig_bounds[1][i])
            lower = min(orig_bounds[0][i], orig_bounds[1][i])
            bounds[i] = [lower, upper]
        return bounds

    def _mesh_res_in_bounds(self, lower, upper, dim):
        """
        Get the mesh resolution in the supplied boundary.

        :param lower: lower boundary.
        :param upper: upper boundary.
        :param dim: Dimension. 0, 1, or 2 for x, y, or z.
        """
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

    def _split_bounds(self, lower, upper, dim):
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

    def _clear_mesh_in_bounds(self, lower, upper, dim):
        for elt in self.mesh_lines[dim]:
            if elt >= lower and elt <= upper:
                self.mesh_lines[dim].remove(elt)

    def _range_union(self, ranges, start_idx=0):
        ranges = sorted(ranges)
        if len(ranges[start_idx:]) <= 1:
            return ranges

        if ranges[start_idx][1] >= ranges[start_idx + 1][0]:
            ranges.append([ranges[start_idx][0], ranges[start_idx + 1][1]])
            del ranges[start_idx : start_idx + 2]
            return self._range_union(ranges[start_idx:])
        else:
            return self._range_union(ranges[start_idx + 1 :])

    def _consolidate_meshed_ranges(self, dim):
        """
        Order meshed ranges and consolidate contiguous ranges.
        """
        self.ranges_meshed[dim] = self._range_union(self.ranges_meshed[dim])

    def _update_ranges(self, lower, upper, dim):
        """
        :param dim: is the dimension: 0, 1, 2 for x, y, or z.
        """
        self.ranges_meshed[dim].append([lower, upper])
        self._consolidate_meshed_ranges(dim)

    def _smooth_mesh_lines(self, dim):
        """
        Ensure adjacent mesh line separations differ by less than the
        smoothness factor.

        If there's enough room between mesh lines, this function will
        recursively add mesh lines corresponding to the largest
        possible separation in line with self.smooth.  When there's
        not enough room, it moves the position of existing lines to be
        in line with smooth.

        TODO should be refactored since logic in if branches are
        basically identical, but switched.

        :param dim: Dimension where mesh should be smoothed.
        """
        for i, pos in enumerate(self.mesh_lines[dim]):
            if i == 0 or i == len(self.mesh_lines[dim]) - 1:
                continue
            left_spacing = pos - self.mesh_lines[dim][i - 1]
            right_spacing = self.mesh_lines[dim][i + 1] - pos
            if (
                left_spacing > (self.smooth * right_spacing)
                and left_spacing - (self.smooth * right_spacing)
                > self.smallest_res / 10
            ):
                ratio = left_spacing / right_spacing
                if i == len(self.mesh_lines[dim]) - 2:
                    # if there's no further mesh spacings to worry
                    # about, this ensures we'll only move the mesh
                    # line when that will satisfy smoothness
                    outer_spacing = right_spacing + (
                        1 / (self.smooth * (self.smooth + 1))
                    )
                else:
                    outer_spacing = self.mesh_lines[dim][i + 2] - (
                        pos + right_spacing
                    )
                # if this condition satisfied, then we can move the
                # current mesh line without violating smoothness
                # elsewhere. To see how I got this condition, imagine
                # spacings are given by a, b and c in order. Spacings
                # on either side of current position are a and b. s
                # gives smooth factor. Move pos by dx to have a and b
                # satisfy s. It's currently above, so move by just
                # enough to satisfy.
                #
                # (a-dx)/(b+dx) = s
                #
                # simultaneously,
                #
                # (b+dx)/c <= s
                #
                # we need the max a where this works. This occurs at
                #
                # (b+dx)/c = s
                #
                # find a. sagemath tells you that
                #
                # a = cs^2 + cs - b
                if (
                    left_spacing
                    <= outer_spacing * self.smooth * (self.smooth + 1)
                    - right_spacing
                ):
                    # adjustment to make left_spacing = smooth * right_spacing
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
                # mesh separation is too small to add smooth *
                # spacing, so instead add it halfway
                elif ratio <= self.smooth * (self.smooth + 1):
                    insort_left(self.mesh_lines[dim], pos - (left_spacing / 2))
                else:
                    insort_left(
                        self.mesh_lines[dim],
                        pos - (self.smooth * right_spacing),
                    )
                self._smooth_mesh_lines(dim)
            elif (
                right_spacing > self.smooth * left_spacing
                and right_spacing - (self.smooth * left_spacing)
                > self.smallest_res / 10
            ):
                ratio = right_spacing / left_spacing
                if i == 1:
                    outer_spacing = left_spacing + (
                        1 / (self.smooth * (self.smooth + 1))
                    )
                else:
                    outer_spacing = (
                        pos - left_spacing - self.mesh_lines[dim][i - 2]
                    )
                if (
                    right_spacing
                    <= outer_spacing * self.smooth * (self.smooth + 1)
                    - left_spacing
                ):
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
                elif ratio <= self.smooth * (self.smooth + 1):
                    insort_left(
                        self.mesh_lines[dim], pos + (right_spacing / 2)
                    )
                else:
                    insort_left(
                        self.mesh_lines[dim],
                        pos + (self.smooth * left_spacing),
                    )
                self._smooth_mesh_lines(dim)

    def _nearest_divisible_res(self, lower, upper, res):
        """
        Return the nearest resolution to @res that evenly subdivides the
        interval [@lower, @upper].

        This is important because it helps prevent adjacent lines from
        bunching up and unnecessarily increasing the simulation time.
        """
        num_divisions = np.round((upper - lower) / res)
        num_divisions = max(num_divisions, 1)
        return (upper - lower) / num_divisions

    def _gen_mesh_in_bounds(self, lower, upper, res, dim, metal=False):
        """
        Add mesh lines within the provided dimensional boundaries.

        :param lower: Lower dimensional boundary.
        :param upper: Upper dimensional boundary.
        :param res: Desired mesh resolution within the provided
            boundary.
        :param dim: Dimension.  0, 1, or 2.
        :param metal: Set to True if this boundary corresponds to a
            metal structure.
        """
        if lower == upper:
            insort_left(self.mesh_lines[dim], lower)
            insort_left(self.const_meshes[dim], lower)
        else:
            [outer_bounds, inner_bounds] = self._split_bounds(
                lower, upper, dim
            )
            for obound in outer_bounds:
                if obound[1] - obound[0] < self.min_lines * res:
                    res = (obound[1] - obound[0]) / self.min_lines
                else:
                    res = self._nearest_divisible_res(
                        obound[0], obound[1], res
                    )
                self.smallest_res = min(self.smallest_res, res)
                j = obound[0]
                while j <= obound[1]:
                    insort_left(self.mesh_lines[dim], j)
                    j += res
                self._update_ranges(obound[0], obound[1], dim)
                if metal:
                    insort_left(self.metal_bounds[dim], obound[0])
                    insort_left(self.metal_bounds[dim], obound[1])
            for ibound in inner_bounds:
                if upper - lower < self.min_lines * res:
                    res = (upper - lower) / self.min_lines
                else:
                    res = self._nearest_divisible_res(
                        ibound[0], ibound[1], res
                    )
                self.smallest_res = min(self.smallest_res, res)
                # only redo the mesh if the desired one is finer than
                # the existing one
                cur_mesh_res = self._mesh_res_in_bounds(
                    ibound[0], ibound[1], dim
                )
                if cur_mesh_res > res and abs(cur_mesh_res - res) > res / 10:
                    self._clear_mesh_in_bounds(ibound[0], ibound[1], dim)
                    j = ibound[0]
                    while j <= ibound[1]:
                        insort_left(self.mesh_lines[dim], j)
                        j += res
                    self._update_ranges(ibound[0], ibound[1], dim)
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
        :param layers: number of conductive layers
        :param sub_epsr: substrate dielectric constant.  dictionary of
            frequency (Hz) and associated dielectric.
        :param sub_rho: volume resistivity (ohm*mm)
        :param layer_sep: separations (in mm) between adjacent copper
            layers.  A list where the first value is the separation
            between the top layer and second layer, etc.  This is
            equivalently the substrate thickness.
        :param layer_thickness: thickness of each conductive layer (in
            mm).  Again proceeds from top to bottom layer.
        """
        self.layers = layers
        self.sub_epsr = sorted(sub_epsr, key=lambda epsr: epsr[0])
        self.sub_rho = sub_rho
        self.sub_kappa = 1 / sub_rho
        self.layer_sep = layer_sep
        self.layer_thickness = layer_thickness

    def epsr_at_freq(self, freq: float):
        """
        Approximate the dielectric at a given frequency given the
        provided epsr values.

        :param freq: frequency of interest (Hz)

        :param returns: dielectric constant
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
    Calculate the microstrip characteristic impedance for a given
    width using Wheeler's equation.  Wheeler's equation can be found
    at:

    https://en.wikipedia.org/wiki/Microstrip#Characteristic_impedance

    :param w: microstrip trace width (mm)
    :param t: trace thickness (mm)
    :param er: substrate relative permittivity
    :param h: substrate height (thickness) (mm)

    :returns: characteristic impedance
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

    :param z0: characteristic impedance (ohm)
    :param t: trace thickness (mm)
    :param er: substrate relative permittivity
    :param h: substrate height (thickness) (mm)
    :param tol: acceptable impedance tolerance (ohm)
    :param guess: an initial guess for the width (mm).  This can
        improve convergence time when the approximate width is known.

    :returns: trace width (mm)
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


class Probe:
    """
    Wrapper class for openems probes.  The main additional feature is
    this class holds onto data.
    """

    def __init__(
        self, name: str, csx: csxcad.ContinuousStructure(), box, p_type=0,
    ):
        """
        :param name: probe name, which becomes filename where data is
            stored.
        :param csx: csx structure to add probe to.
        :param box: 2D list, outer dim=2, inner dim=3.  Outer
            dimensions is start and stop of box, and inner dimension
            are the x, y, and z coordinates, respectively.
        :param p_type: 0 if voltage probe, 1 if current probe.
        """
        self.data = None
        self.f_data = None
        self.box = box
        self.csx = csx
        self.name = name
        self.p_type = p_type
        self.csxProbe = self.csx.AddProbe(self.name, self.p_type)
        self.csxProbe.AddBox(start=self.box[0], stop=self.box[1])

    def calc_frequency_data(self, sim_path, freq, signal_type="pulse"):
        self.data = UI_data([self.name], sim_path, freq, signal_type)
        self.f_data = self.data.ui_f_val[0]


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
        fcsx: str = None,
        fvtr_dir: str = None,
        efield: bool = False,
    ):
        """
        :param pcb: PCB object
        :param f0: center frequency (Hz)
        :param fc: corner frequency, given as difference from f0 (Hz)
        :param z0_ref: desired impedance (ohm)
        :param microstrip_width: microstrip trace width (mm).  If
            ommitted, an analytical best guess for z0_ref and f0 will
            be used (see `wheeler_z0_width`).
        :param microstrip_len: microstrip trace length (mm).  This
            value isn't critical but does need to be large enough for
            the openems simulation to work.  If unsure, stick with the
            default.
        :param substrate_width: substrate width (mm)
        :param substrate_len: substrate length (mm).  Must be at least
            as large as the microstrip length.  If unsure, use the
            default.
        :param fcsx: CSX file name to write CSXCAD structure to.  If
            ommitted, do not save the file.  Omitting this does not
            prevent you from viewing the CSX structure, just from
            saving it for later.  Relative to current directory.
        :param fvtr_dir: directory for VTR E-field dump files.  Files
            will be prefixed with 'Et_'.  If ommitted, E-field can
            still be viewed unless `efield` is set to false in which
            case this parameter has no effect.
        :param efield: dump E-field time values for viewing.
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
        if fcsx is None:
            self.fcsx = tempfile.mkstemp()[1]
        else:
            self.fcsx = os.path.abspath(fcsx)

        self.efield = efield
        if fvtr_dir is None:
            if efield:
                self.fvtr_dir = tempfile.mkdtemp()
        else:
            self.fvtr_dir = os.path.abspath(fvtr_dir)
            if os.path.exists(self.fvtr_dir):
                rmtree(self.fvtr_dir)
            os.mkdir(self.fvtr_dir)

        self.mesh = None
        self.csx = None
        self.vprobes = [None]
        self.iprobes = [None]
        # track simulation steps so that we don't reperform steps
        self.sim_done = False
        self.csx_done = False
        self.freq_bins = None
        self.fdtd = None

    def sim(
        self, num_freq_bins: int = 501, zero_trace_height: bool = False,
    ):
        """
        Run an openems simulation for the current microstrip structure
        and write the results to a file.

        :param num_freq_bins: number of frequency bins.  More
            frequency bins means longer simulation time but possibly
            greater accuracy.
        :param zero_trace_height: approximate trace height as 0.  This
            can greatly reduce simulation time but may reduce
            accuracy.
        """
        self.gen_csx(zero_trace_height=zero_trace_height)
        tmpdir = tempfile.mkdtemp()
        self.fdtd.Run(tmpdir, cleanup=True)

        self.freq_bins = np.linspace(
            self.f0 - self.fc, self.f0 + self.fc, num_freq_bins
        )

        for i, _ in enumerate(self.vprobes):
            self.vprobes[i].calc_frequency_data(
                sim_path=tmpdir, freq=self.freq_bins
            )

        for i, _ in enumerate(self.iprobes):
            self.iprobes[i].calc_frequency_data(
                sim_path=tmpdir, freq=self.freq_bins
            )

        self.calc_params()
        self.sim_done = True
        return np.absolute(self.Z_ref)

    def gen_csx(self, zero_trace_height: bool = False) -> None:
        """
        Generate CSX structure.

        :param zero_trace_height: approximate trace height as 0.  This
            can greatly reduce simulation time but may reduce
            accuracy.  Args: num_probes: zero_trace_height:
        """
        # dimensions and parameters
        unit = 1e-3
        lmin = C0 / ((self.f0 + self.fc) * unit)
        if zero_trace_height:
            trace_height = 0
        else:
            trace_height = self.pcb.layer_thickness[0]

        # structures
        self.csx = csxcad.ContinuousStructure()
        microstrip = self.csx.AddMetal("PEC")
        microstrip.AddBox(
            priority=10,
            start=[-self.microstrip_len / 2, -self.microstrip_width / 2, 0],
            stop=[
                self.microstrip_len / 2,
                self.microstrip_width / 2,
                trace_height,
            ],
        )
        substrate = self.csx.AddMaterial(
            "substrate",
            epsilon=self.pcb.epsr_at_freq(self.f0),
            kappa=self.pcb.sub_kappa,
        )
        substrate.AddBox(
            priority=0,
            start=[
                -self.microstrip_len / 2,
                -self.substrate_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[self.microstrip_len / 2, self.substrate_width / 2, 0],
        )

        # simulation
        self.fdtd = openems.openEMS(EndCriteria=1e-5)
        self.fdtd.SetCSX(self.csx)
        self.fdtd.SetGaussExcite(self.f0, self.fc)
        self.fdtd.SetBoundaryCond(
            ["PML_8", "PML_8", "MUR", "MUR", "PEC", "MUR"]
        )

        self.mesh = Mesh(
            self.csx,
            lmin,
            mres=1 / 20,
            sres=1 / 10,
            smooth=1.4,
            unit=unit,
            min_lines=9,
            expand_bounds=[0, 0, 10, 10, 0, 10],
        )
        self.mesh.generate_mesh()

        _, port_xpos = self.mesh.nearest_mesh_line(0, -self.microstrip_len / 4)
        self.fdtd.AddLumpedPort(
            port_nr=0,
            R=self.z0_ref,
            start=[
                port_xpos,
                -self.microstrip_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[port_xpos, self.microstrip_width / 2, 0],
            p_dir="z",
            excite=1,
            priority=999,
        )

        self.gen_probes(trace_height=trace_height)

        # E-field recording
        if self.efield:
            Et = self.csx.AddDump(
                os.path.join(self.fvtr_dir, "Et_"), file_type=0
            )
            start = [
                -self.microstrip_len / 2,
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
        self.csx.Write2XML(self.fcsx)
        self.csx_done = True

    def gen_probes(self, trace_height):
        mid_idx, mid_xpos = self.mesh.nearest_mesh_line(0, 0)
        vprobe_x_pos = [
            self.mesh.mesh_lines[0][mid_idx - 1],
            mid_xpos,
            self.mesh.mesh_lines[0][mid_idx + 1],
        ]
        iprobe_x_pos = [
            (vprobe_x_pos[0] + vprobe_x_pos[1]) / 2,
            (vprobe_x_pos[1] + vprobe_x_pos[2]) / 2,
        ]

        self.vprobes = [None for i in range(len(vprobe_x_pos))]
        self.iprobes = [None for i in range(len(iprobe_x_pos))]

        for i, _ in enumerate(self.vprobes):
            self.vprobes[i] = Probe(
                name="ut_" + str(i),
                csx=self.csx,
                box=[
                    [vprobe_x_pos[i], 0, -self.pcb.layer_sep[0]],
                    [vprobe_x_pos[i], 0, 0],
                ],
                p_type=0,
            )

        for i, _ in enumerate(self.iprobes):
            self.iprobes[i] = Probe(
                name="it_" + str(i),
                csx=self.csx,
                box=[
                    [iprobe_x_pos[i], -self.microstrip_width / 2, 0],
                    [
                        iprobe_x_pos[i],
                        self.microstrip_width / 2,
                        trace_height,
                    ],
                ],
                p_type=1,
            )

    def calc_params(self):
        Et = self.vprobes[1].f_data
        Ht = 0.5 * (self.iprobes[0].f_data + self.iprobes[1].f_data)
        unit = 1e-3
        dEt = (self.vprobes[2].f_data - self.vprobes[0].f_data) / (
            unit * (self.vprobes[2].box[0][0] - self.vprobes[0].box[0][0])
        )
        dHt = (self.iprobes[1].f_data - self.iprobes[0].f_data) / (
            unit * (self.iprobes[1].box[0][0] - self.iprobes[0].box[0][0])
        )
        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[
            np.real(beta) < 0
        ] *= -1  # determine correct sign (unlike the paper)
        self.beta = beta

        # determine ZL
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))

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


def microstrip_sweep_width(
    pcb: PCB,
    f0: float,
    fc: float,
    z0_ref: float,
    width: float = None,
    width_dev_factor: float = 0.1,
    num_points: int = 11,
) -> List[Tuple[float, float]]:
    """
    Calculate microstrip characteristic impedance as a function of
    trace width.

    :param pcb: PCB object
    :param f0: center frequency (Hz)
    :param fc: corner frequency, given as difference from f0 (Hz)
    :param z0_ref: desired impedance (ohm)
    :param width: center width (mm).  If ommitted, find the best guess
        width for f0 analytically.
    :param width_dev_factor: determines width sweep bounds
        [width*(1-wdf),width*(1+wdf)]
    :param num_points: number of width values to calculate.  This
        value has a significant effect on simulation time since each
        point is calculated in its own thread.  Therefore, to minimize
        computation time, its recommended to choose some multiple of
        the number of cores available on the simulation machine.  The
        total simulation time will be roughly equal to the time it
        takes to compute 1 point times the ratio of the number of
        points to number of machine cores.
    :param fout: filename where results should be written.  relative
        to current dir.

    :returns: A list of tuples, where the first tuple element is a
        width and the second is the corresponding impedance value.
    """
    if width is None:
        width = wheeler_z0_width(
            z0=z0_ref,
            t=pcb.layer_thickness[0],
            er=pcb.epsr_at_freq(f0),
            h=pcb.layer_sep[0],
        )

    widths = np.linspace(
        width * (1 - width_dev_factor),
        width * (1 + width_dev_factor),
        num_points,
    )
    microstrips = [None for i in range(num_points)]
    for i, width in enumerate(widths):
        microstrips[i] = Microstrip(
            pcb=pcb, f0=f0, fc=fc, z0_ref=z0_ref, microstrip_width=width
        )

    pool = Pool(nodes=11)
    freq_bins = 501
    func = partial(
        Microstrip.sim, num_freq_bins=freq_bins, zero_trace_height=False,
    )
    z0s = [None for i in range(num_points)]
    z0s = list(pool.map(func, microstrips))

    ret_vals = [[None, None] for i in range(num_points)]
    for i, _ in enumerate(z0s):
        f0_bin_idx = int(freq_bins / 2)
        z0 = z0s[i]
        ret_vals[i][0] = widths[i]
        ret_vals[i][1] = z0[f0_bin_idx]

    return ret_vals
