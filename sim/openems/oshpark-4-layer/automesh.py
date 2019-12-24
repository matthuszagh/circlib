from bisect import insort_left, bisect_left, bisect_right
import numpy as np


class AutoMesh:
    def __init__(
        self,
        csx,
        lmin,
        simbox,
        mres=1 / 20,
        sres=1 / 10,
        smooth=1.5,
        unit=1e-3,
        min_lines=10,
    ):
        """
        @csx is the CSXCAD structure (return value of
        CSXCAD.ContinuousStructure()).
        @lmin is the minimum wavelength associated with the expected
        frequency.
        @simbox lists the simulation box parameters in the form
        [[xmin, xmax], [ymin, ymax], [zmin, zmax]].
        @mres is the metal resolution, specified as a factor of @lmin.
        @mres is the substrate resolution, specified as a factor of @lmin.
        @smooth is the factor by which adjacent cells are allowed to differ in
        size.
        @unit is the mesh size unit, which defaults to mm.
        @min_lines is the minimum number of mesh lines for a primitive's
        dimensional length, unless the length of that primitive's dimension
        is precisely 0.
        """
        self.csx = csx
        self.lmin = lmin
        self.simbox = simbox
        self.mres = mres * self.lmin
        self.sres = sres * self.lmin
        self.smooth = smooth
        self.min_lines = min_lines
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
        exceptions to this: (1) if a dimension of a meterial has length 0
        (e.g. a planar metal sheet) a single mesh line is placed exactly
        on that line, and (2) a nonzero length material must have a
        minimum of 10 mesh lines (TODO this choice should be
        evaluated). The 2nd exception will not violate the third's rule or
        smoothness (that adjacent mesh lines not differ in separation by
        more than a factor of 1.5, TODO which should also be evaluated).
        """
        # add primitive meshes
        for prim in self.prims:
            bounds = self.GetPrimBounds(prim)
            if self.PrimMetalp(prim):
                for i in range(3):
                    self.GenMeshInBounds(
                        bounds[i][0], bounds[i][1], self.mres, i, metal=True
                    )
            else:
                for i in range(3):
                    self.GenMeshInBounds(
                        bounds[i][0], bounds[i][1], self.sres, i, metal=False
                    )

        # add simulation box mesh
        for i in range(3):
            self.GenMeshInBounds(
                self.simbox[i][0], self.simbox[i][1], self.sres, i, metal=False
            )
        # remove unintended, tightly spaced meshes
        for dim in range(3):
            self.RemoveTightMeshLines(dim)

        # enforce thirds rule
        for dim in range(3):
            self.EnforceThirds(dim)

        # smooth mesh
        for dim in range(3):
            self.SmoothMeshLines(dim)
        # set calculated mesh lines
        self.AddAllMeshLines()

    # TODO should ensure that inserted mesh lines are not at metal boundaries
    def EnforceThirds(self, dim):
        for i, pos in enumerate(self.mesh_lines[dim]):
            if pos in self.metal_bounds[dim]:
                # at lower boundary
                if i == 0:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos + (self.mres / 3))
                    self.EnforceThirds(dim)
                # at upper boundary
                elif i == len(self.mesh_lines[dim]) - 1:
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos - (self.mres / 3))
                    self.EnforceThirds(dim)
                else:
                    spacing_left = pos - self.mesh_lines[dim][i - 1]
                    spacing_right = self.mesh_lines[dim][i + 1] - pos
                    del self.mesh_lines[dim][i]
                    # metal-metal boundary
                    if spacing_left == spacing_right:
                        insort_left(
                            self.mesh_lines[dim], pos - (self.mres / 3)
                        )
                        insort_left(
                            self.mesh_lines[dim], pos + (self.mres / 3)
                        )
                    elif spacing_left < spacing_right:
                        insort_left(
                            self.mesh_lines[dim], pos - (self.mres / 3)
                        )
                        insort_left(
                            self.mesh_lines[dim], pos + (2 * self.mres / 3)
                        )
                    else:
                        insort_left(
                            self.mesh_lines[dim], pos - (2 * self.mres / 3)
                        )
                        insort_left(
                            self.mesh_lines[dim], pos + (self.mres / 3)
                        )
                    self.EnforceThirds(dim)

    def RemoveTightMeshLines(self, dim):
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
                self.RemoveTightMeshLines(dim)
            # we have to check whether these are zero-dimension
            # structures before deleting them.
            elif (
                pos - last_pos < self.smallest_res
                and pos not in self.const_meshes[dim]
                and last_pos not in self.const_meshes[dim]
            ):
                if last_pos not in self.const_meshes[dim]:
                    del self.mesh_lines[dim][i - 1]
                else:
                    del self.mesh_lines[dim][i]
                self.RemoveTightMeshLines(dim)
            else:
                last_pos = pos

    def AddAllMeshLines(self):
        for dim in range(3):
            for line in self.mesh_lines[dim]:
                self.mesh.AddLine(dim, line)

    def GetMesh(self):
        return self.mesh

    def PrimMetalp(self, prim):
        return prim.GetProperty().GetTypeString() == "Metal"

    def GetPrimBounds(self, prim):
        orig_bounds = prim.GetBoundBox()
        bounds = [[None, None], [None, None], [None, None]]
        for i in range(3):
            upper = max(orig_bounds[0][i], orig_bounds[1][i])
            lower = min(orig_bounds[0][i], orig_bounds[1][i])
            bounds[i] = [lower, upper]
        return bounds

    def MeshResInBounds(self, lower, upper, dim):
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

    def SplitBounds(self, lower, upper, dim):
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
                outin_ranges[1].append([lower, upper_mesh])
                if upper > upper_mesh:
                    lower = upper_mesh
                    continue
                else:
                    return outin_ranges
        if lower < upper:
            outin_ranges[0].append([lower, upper])

        return outin_ranges

    def ClearMeshInBounds(self, lower, upper, dim):
        for elt in self.mesh_lines[dim]:
            if elt >= lower and elt <= upper:
                self.mesh_lines.remove(elt)

    def RangeUnion(self, ranges, start_idx=0):
        ranges = sorted(ranges)
        if len(ranges[start_idx:]) <= 1:
            return ranges

        if ranges[start_idx][1] >= ranges[start_idx + 1][0]:
            ranges.append([ranges[start_idx][0], ranges[start_idx + 1][1]])
            del ranges[start_idx : start_idx + 2]
            return self.RangeUnion(ranges[start_idx:])
        else:
            return self.RangeUnion(ranges[start_idx + 1 :])

    def ConsolidateMeshedRanges(self, dim):
        """
        Order meshed ranges and consolidate contiguous ranges.
        """
        self.ranges_meshed[dim] = self.RangeUnion(self.ranges_meshed[dim])

    def UpdateRanges(self, lower, upper, dim):
        """
        @dim is the dimension: 0, 1, 2 for x, y, or z.
        """
        self.ranges_meshed[dim].append([lower, upper])
        self.ConsolidateMeshedRanges(dim)

    def SmoothMeshLines(self, dim):
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
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos - adj)
                else:
                    insort_left(self.mesh_lines[dim], pos - (left_spacing / 2))
                self.SmoothMeshLines(dim)
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
                    del self.mesh_lines[dim][i]
                    insort_left(self.mesh_lines[dim], pos + adj)
                else:
                    insort_left(
                        self.mesh_lines[dim], pos + (right_spacing / 2)
                    )
                self.SmoothMeshLines(dim)

    def NearestDivisibleRes(self, lower, upper, res):
        """
        Return the nearest resolution to @res that evenly subdivides the
        interval [@lower, @upper].

        This is important because it helps prevent adjacent lines from
        bunching up and unnecessarily increasing the simulation time.
        """
        num_divisions = np.round((upper - lower) / res)
        num_divisions = max(num_divisions, 1)
        return (upper - lower) / num_divisions

    def GenMeshInBounds(self, lower, upper, res, dim, metal=False):
        if lower == upper:
            insort_left(self.mesh_lines[dim], lower)
            insort_left(self.const_meshes[dim], lower)
        else:
            [outer_bounds, inner_bounds] = self.SplitBounds(lower, upper, dim)
            for obound in outer_bounds:
                if upper - lower < self.min_lines * res:
                    res = (upper - lower) / self.min_lines
                else:
                    res = self.NearestDivisibleRes(obound[0], obound[1], res)
                self.smallest_res = min(self.smallest_res, res)
                j = obound[0]
                while j <= obound[1]:
                    insort_left(self.mesh_lines[dim], j)
                    j += res
                self.UpdateRanges(obound[0], obound[1], dim)
                if metal:
                    insort_left(self.metal_bounds[dim], obound[0])
                    insort_left(self.metal_bounds[dim], obound[1])
            for ibound in inner_bounds:
                if upper - lower < self.min_lines * res:
                    res = (upper - lower) / self.min_lines
                else:
                    res = self.NearestDivisibleRes(ibound[0], ibound[1], res)
                self.smallest_res = min(self.smallest_res, res)
                # only redo the mesh if the desired one is finer than
                # the existing one
                cur_mesh_res = self.MeshResInBounds(ibound[0], ibound[1], dim)
                if cur_mesh_res > res and abs(cur_mesh_res - res) > res / 10:
                    self.ClearMeshInBounds(ibound[0], ibound[1], dim)
                    j = ibound[0]
                    while j <= ibound[1]:
                        insort_left(self.mesh_lines[dim], j)
                        j += res
                    self.UpdateRanges(ibound[0], ibound[1], dim)
                if metal:
                    insort_left(self.metal_bounds[dim], ibound[0])
                    insort_left(self.metal_bounds[dim], ibound[1])
