from concurrent.futures import ThreadPoolExecutor

import numpy as np

from LoopStructural.utils.maths import gradient_from_tetrahedron, regular_tetraherdron_for_points

from ....modelling.features import FeatureType, StructuralFrame
from ....modelling.features.fault._fault_function import BaseFault, BaseFault3D, FaultDisplacement
from ....modelling.features.fault._fault_function_feature import (
    FaultDisplacementFeature,
)
from ....utils import NegativeRegion, PositiveRegion, getLogger
from .._evaluation_cache import bump_revision, values_differ

logger = getLogger(__name__)

use_threads = False  # there is a bug with threading causing the interpolation dictionary to change size during iterations
# for now 19/4/24 LG is turning threading off. This will slow down evaluation of models a bit.


class FaultSegment(StructuralFrame):
    """
    Class for representing a slip event of a fault
    """

    def __init__(
        self,
        name: str,
        features: list,
        faultfunction=None,
        steps=10,
        displacement=1.0,
        fold=None,
        model=None,
    ):
        """
        A slip event of a fault

        Parameters
        ----------
        faultframe : FaultFrame
            the fault frame defining the faut geometry
        faultfunction : function/lambda function
            optional displacement function for spatially variable fault displacement
        steps : int
            how many integration steps for faults
        kwargs
        """
        StructuralFrame.__init__(self, name=name, features=features, fold=fold, model=model)
        self.type = FeatureType.FAULT
        self.displacement = displacement
        self._faultfunction = BaseFault().fault_displacement
        self.steps = steps
        self.regions = []
        self.faults_enabled = True

        self.builder = None
        self.splay = {}
        self.abut = {}
        self.fault_offset = 0.0

    @property
    def displacement(self):
        return self._displacement

    @displacement.setter
    def displacement(self, value):
        changed = values_differ(getattr(self, "_displacement", None), value)
        self._displacement = value
        if changed:
            bump_revision(self)

    @property
    def steps(self):
        return self._steps

    @steps.setter
    def steps(self, value):
        changed = values_differ(getattr(self, "_steps", None), value)
        self._steps = value
        if changed:
            bump_revision(self)

    @property
    def faultfunction(self):
        return self._faultfunction

    @faultfunction.setter
    def faultfunction(self, value):
        try:
            if callable(value):
                self._faultfunction = value
            if issubclass(FaultDisplacement, type(value)):
                self._faultfunction = value

            elif isinstance(value, str) and value == "BaseFault":
                self._faultfunction = BaseFault().fault_displacement
            elif isinstance(value, str) and value == "BaseFault3D":
                self._faultfunction = BaseFault3D().fault_displacement
            else:
                raise ValueError("Fault function must be a function or BaseFault")
        finally:
            # The function can change before an error. The revision must grow also then.
            bump_revision(self)

    @property
    def fault_normal_vector(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_normal_vector

    @property
    def fault_slip_vector(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_slip_vector

    @property
    def fault_strike_vector(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_strike_vector

    @property
    def fault_minor_axis(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_minor_axis

    @property
    def fault_major_axis(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_major_axis
    
    @property
    def fault_intermediate_axis(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_intermediate_axis

    @property
    def fault_centre(self):
        if self.builder is None:
            raise ValueError("Fault builder not set")
        return self.builder.fault_centre

    @property
    def displacementfeature(self):
        return FaultDisplacementFeature(self, self.faultfunction, name=self.name, model=self.model)

    def fault_ellipsoid(self, **kwargs):
        try:
            import pyvista as pv

            if self.model is None:
                pts = self.fault_centre[None, :]
            else:
                pts = self.model.rescale(self.fault_centre[None, :], inplace=False)
            # pts = self.fault_centre[None, :]
            fault_ellipsoid = pv.PolyData(pts)
            # fault_ellipsoid = pv.PolyData(
            #     self.model.rescale(self.fault_centre[None, :], inplace=False)
            # )
            fault_ellipsoid["norm"] = self.fault_normal_vector[None, :]
            fault_ellipsoid['norm'] /= np.linalg.norm(fault_ellipsoid['norm'], axis=1)[:, None]
            geom = pv.ParametricEllipsoid(
                self.fault_minor_axis,
                self.fault_major_axis,
                self.fault_intermediate_axis,
            )
            ellipsoid = fault_ellipsoid.glyph(geom=geom, **kwargs)
            return ellipsoid
        except ImportError:
            logger.error("pyvista not installed")
            return None

    def set_fault_offset(self, offset: float):
        self.fault_offset = offset

    def set_model(self, model):
        """
        Link a geological model to the feature

        Parameters
        ----------
        model - GeologicalModel

        Returns
        -------

        """
        self.model = model

    def set_displacement(self, displacement, scale=True):
        """
        Set the fault displacement to a new value

        Parameters
        ----------
        displacement - double
        scale - boolean

        Returns
        -------

        """
        if scale and self.model is not None:
            self.displacement = displacement / self.model.scale_factor
            bump_revision(self)
        elif not scale:
            self.displacement = displacement
            bump_revision(self)
        else:
            logger.warning("Displacement not updated")

    def evaluate(self, locations):
        """
        Evaluate which side of fault

        Parameters
        ----------
        locations numpy array
            location to evaluate

        Returns
        -------
            boolean array true if on hanging wall, false if on footwall

        """
        v = self.__getitem__(0).evaluate_value(locations)
        v[~np.isnan(v)] = v[~np.isnan(v)] > 0
        v[np.isnan(v)] = 0
        return v.astype(bool)

    def inside_volume(self, locations, threshold=0.001):
        # v = self.faultframe.evaluate_value(locations)
        v = self.evaluate_displacement(locations) / self.displacement
        v[np.isnan(v)] = 0
        return np.abs(v) > threshold
        # return np.all(np.logical_and(v > -1,v<1),axis=1)

    def evaluate_value(self, locations, ignore_regions=False):
        """
        Return the value of the fault surface scalar field

        Parameters
        ----------
        locations - numpy array
            location to evaluate scalar field

        Returns
        -------

        """
        # v = np.zeros(locations.shape[0])
        # v[:] = np.nan
        # mask =
        # mask = np.zeros(locations.shape[0]).astype(bool)
        # mask[:] = True
        # # check regions
        # for r in self.regions:
        #     try:
        #         mask = np.logical_and(mask, r(locations))
        #     except:
        #         logger.error("nan slicing")
        # v[mask] = self.__getitem__(0).evaluate_value(locations[mask, :])
        return self.__getitem__(0).evaluate_value(locations, ignore_regions=ignore_regions)

    def ellipsoid(self):
        pass

    def mean(self):
        return self.__getitem__(0).mean()

    def max(self):
        return self.__getitem__(0).max()

    def min(self):
        return self.__getitem__(0).min()

    def evaluate_gradient(self, locations):
        """
        Return the fault slip direction at the location

        Parameters
        ----------
        locations - numpy array Nx3


        Returns
        -------

        """
        v = np.zeros(locations.shape)
        v[:, :] = np.nan
        mask = np.zeros(locations.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            try:
                mask = np.logical_and(mask, r(locations))
            except (ValueError, IndexError) as e:
                logger.error(f"nan slicing: {e}")
        # need to scale with fault displacement
        v[mask, :] = self.__getitem__(1).evaluate_gradient(locations[mask, :])
        v[mask, :] /= np.linalg.norm(v[mask, :], axis=1)[:, None]
        scale = self.displacementfeature.evaluate_value(locations[mask, :])
        v[mask, :] *= scale[:, None] * self.displacement
        return v

    def evaluate_displacement(self, points):
        newp = np.copy(points).astype(float)
        # evaluate fault function for all
        # points then define mask for only points affected by fault
        gx = None
        gy = None
        gz = None
        if use_threads:
            with ThreadPoolExecutor(max_workers=8) as executor:
                # all of these operations should be independent so
                # just run as different threads
                gx_future = executor.submit(self.__getitem__(0).evaluate_value, newp)
                gy_future = executor.submit(self.__getitem__(1).evaluate_value, newp)
                gz_future = executor.submit(self.__getitem__(2).evaluate_value, newp)
                gx = gx_future.result()
                gy = gy_future.result()
                gz = gz_future.result()
        else:
            gx = self.__getitem__(0).evaluate_value(newp)
            gy = self.__getitem__(1).evaluate_value(newp)
            gz = self.__getitem__(2).evaluate_value(newp)
        d = np.zeros(gx.shape)
        mask = np.logical_and(~np.isnan(gx), ~np.isnan(gy))
        mask = np.logical_and(mask, ~np.isnan(gz))
        d[~mask] = 0
        gx_mask = np.zeros_like(mask, dtype=bool)
        gx_mask[mask] = gx[mask] > 0
        d[gx_mask] = 1.0
        if self.faultfunction is not None:
            d[mask] = self.faultfunction(gx[mask] + self.fault_offset, gy[mask], gz[mask])
        return d * self.displacement

    def _shared_frame_members(self):
        """Get the frame members that can share one grid lookup.

        The fast path is only safe when all of these are true for the three
        members of the frame:

        - each member is a plain GeologicalFeature with a discrete interpolator
          on a StructuredGrid that has the uniform cell lookup,
        - the grids have equal geometry and the interpolators have equal
          transforms,
        - no member has faults,
        - the regions of the members are the same objects, so one mask is
          valid for all members.

        Returns
        -------
        list or None
            the three members, or None if the legacy calls are needed
        """
        from loop_common.supports import StructuredGrid
        from loop_interpolation import DiscreteInterpolator

        from .._geological_feature import GeologicalFeature

        try:
            members = [self.__getitem__(i) for i in range(3)]
        except (IndexError, KeyError, TypeError):
            return None
        reference = None
        for member in members:
            if type(member).evaluate_value is not GeologicalFeature.evaluate_value:
                return None
            if type(member).evaluate_gradient is not GeologicalFeature.evaluate_gradient:
                return None
            if len(member.faults) > 0:
                return None
            interpolator = member.interpolator
            if not isinstance(interpolator, DiscreteInterpolator):
                return None
            support = interpolator.support
            if not isinstance(support, StructuredGrid) or not support._has_uniform_lookup():
                return None
            box = interpolator.bounding_box
            description = (
                np.asarray(support.origin),
                np.asarray(support.nsteps),
                np.asarray(support.step_vector),
                None if box is None else np.asarray(box.local_origin),
                None if box is None else np.asarray(box.local_rotation),
                member.regions,
            )
            if reference is None:
                reference = description
                continue
            for index, (a, b) in enumerate(zip(reference, description)):
                if index == 5:
                    if len(a) != len(b) or any(x is not y for x, y in zip(a, b)):
                        return None
                elif a is None or b is None:
                    if a is not b:
                        return None
                elif not np.array_equal(a, b):
                    return None
        return members

    def _evaluate_frame(self, points, need_gradient=False, need_values=True, need_gx=True):
        """Evaluate the fault frame at the points.

        If the three members have equal StructuredGrid supports, one grid
        lookup gives the values of all members and the gradient of coordinate
        1. Other frames use the calls of the members (the legacy path).
        The lookup is not kept after this call.

        Parameters
        ----------
        points : np.ndarray
            Nx3 array of points
        need_gradient : bool
            add the gradient of coordinate 1 to the result. This gradient is not
            normalised and it ignores the regions.
        need_values : bool
            if false, the values are not evaluated and are None
        need_gx : bool
            if false, the values of coordinate 0 are not evaluated and are None

        Returns
        -------
        tuple
            (gx, gy, gz) or (gx, gy, gz, gradient)
        """
        points = np.asarray(points)
        members = self._shared_frame_members() if len(points) > 0 else None
        if members is None:
            return self._evaluate_frame_legacy(points, need_gradient, need_values, need_gx)
        for member in members:
            if member.builder is not None:
                member.builder.up_to_date()
        interpolators = [member.interpolator for member in members]
        for interpolator in interpolators:
            interpolator.update()
        support = interpolators[0].support
        local = interpolators[0]._project_points(points)
        lookup = support._lookup_points(local)
        values = [None, None, None]
        if need_values:
            mask = members[0]._calculate_mask(points)
            nan_rows = np.any(np.isnan(np.asarray(local)), axis=1)
            for i in range(3):
                if i == 0 and not need_gx:
                    continue
                v = support._evaluate_value_from_lookup(lookup, interpolators[i].c)
                # The interpolator gives 0 for points that are not finite.
                v[nan_rows] = 0.0
                v[~mask] = np.nan
                values[i] = v
        if not need_gradient:
            return tuple(values)
        gradient = support._evaluate_gradient_from_lookup(lookup, interpolators[1].c)
        gradient = interpolators[1]._reproject_vectors(gradient)
        return (*values, gradient)

    def _evaluate_frame_legacy(self, points, need_gradient, need_values, need_gx):
        """Evaluate the frame with one call for each member and gradient."""
        member0, member1, member2 = (self.__getitem__(i) for i in range(3))
        gx = gy = gz = g = None
        if use_threads:
            with ThreadPoolExecutor(max_workers=8) as executor:
                # all of these operations should be
                # independent so just run as different threads
                futures = {}
                if need_values:
                    if need_gx:
                        futures["gx"] = executor.submit(member0.evaluate_value, points)
                    futures["gy"] = executor.submit(member1.evaluate_value, points)
                    futures["gz"] = executor.submit(member2.evaluate_value, points)
                if need_gradient:
                    futures["g"] = executor.submit(
                        member1.evaluate_gradient, points, ignore_regions=True
                    )
                gx = futures["gx"].result() if "gx" in futures else None
                gy = futures["gy"].result() if "gy" in futures else None
                gz = futures["gz"].result() if "gz" in futures else None
                g = futures["g"].result() if "g" in futures else None
        else:
            if need_values:
                if need_gx:
                    gx = member0.evaluate_value(points)
                gy = member1.evaluate_value(points)
                gz = member2.evaluate_value(points)
            if need_gradient:
                g = member1.evaluate_gradient(points, ignore_regions=True)
        if need_gradient:
            return gx, gy, gz, g
        return gx, gy, gz

    def apply_to_points(self, points, reverse=False, precomputed_gx=None):
        """
        Unfault the array of points

        Parameters
        ----------
        points - numpy array Nx3
        precomputed_gx - optional pre-evaluated gx values (same points, avoids duplicate evaluation)

        Returns
        -------

        """
        logger.info(f'Applying fault {self.name} to points {points.shape}')
        steps = self.steps
        newp = np.copy(points).astype(float)
        # evaluate fault function for all points
        # then define mask for only points affected by fault
        # gx may be supplied by caller to avoid re-evaluation (precomputed from region check)
        if precomputed_gx is not None:
            gx = precomputed_gx
            _, gy, gz = self._evaluate_frame(newp, need_gx=False)
        else:
            gx, gy, gz = self._evaluate_frame(newp)
        d = np.zeros(gx.shape)
        mask = np.logical_and(~np.isnan(gx), ~np.isnan(gy))
        mask = np.logical_and(mask, ~np.isnan(gz))
        d[~mask] = 0
        gx_mask = np.zeros_like(mask, dtype=bool)
        gx_mask[mask] = gx[mask] > 0
        d[gx_mask] = 1.0
        if self.faultfunction is not None:
            d[mask] = self.faultfunction(gx[mask], gy[mask], gz[mask])
        mask = np.abs(d) > 0.0

        d *= self.displacement
        if reverse:
            d *= -1.0
        # calculate the fault frame for the evaluation points
        for _i in range(steps):
            if _i == 0:
                # Reuse gx/gy/gz from the initial full-array evaluation above — newp[mask] hasn't
                # moved yet on the first iteration, so values are identical.
                gx_m = gx[mask]
                gy_m = gy[mask]
                gz_m = gz[mask]
                g = self._evaluate_frame(newp[mask, :], need_gradient=True, need_values=False)[3]
            else:
                gx_m, gy_m, gz_m, g = self._evaluate_frame(newp[mask, :], need_gradient=True)
            gx, gy, gz = gx_m, gy_m, gz_m  # alias for block below
            # # get the fault frame val/grad for the points
            # determine displacement magnitude, for constant displacement
            # hanging wall should be > 0
            d = np.zeros(gx.shape)
            mask2 = np.logical_and(~np.isnan(gx), ~np.isnan(gy))
            mask2 = np.logical_and(mask2, ~np.isnan(gz))
            d[~mask2] = 0
            gx_mask2 = np.zeros_like(mask2, dtype=bool)
            gx_mask2[mask2] = gx[mask2] > 0
            # d[~np.isnan(gx)][gx[~np.isnan(gx)]>0] = 1
            d[gx_mask2] = 1.0
            # d[mask2][gx[mask2] < 0] = 0.
            # d[gx < 0] = 0.
            if self.faultfunction is not None:
                d[mask2] = self.faultfunction(gx[mask2], gy[mask2], gz[mask2])
            d *= self.displacement
            # normalise when length is >0
            g_mag = np.zeros(g.shape[0])
            g_mag[mask2] = np.linalg.norm(g[mask2], axis=1)
            # g_mag = np.linalg.norm(g[mask2], axis=1)
            g[g_mag > 0.0] /= g_mag[g_mag > 0, None]
            # multiply displacement vector by the displacement magnitude for
            # step

            g *= (1.0 / steps) * d[:, None]
            # newp[mask, :].copy()
            # apply displacement
            newp[mask, :] += g

        return newp

    def apply_to_vectors(self, vector: np.ndarray, scale_parameter: float = 1.0) -> np.ndarray:
        """Rotates the vector through the fault displacement field

        Parameters
        ----------
        vector : np.ndarray
            Nx6 array of vectors x,y,z,vx,vy,vz
        """

        # define a regular tetrahedron as a mask to apply to every vector observation.
        # the regular tetrahedron is located at the location of the vector datapoint
        # and the corner values of the tetrahedron are defined by setting the barycenter
        # of the tetrahedron as 0 and then using the vector to calculate the values of the
        # corner assuming the tetra is a P1 tetra with a linear gradient.
        # a scaling parameter is used to determine the size of the tetrahedron
        # the nodes of the tetrahedron are then restored by the fault and then gradient is
        # recalculated using the updated node positions but the original corner values

        tetrahedron = regular_tetraherdron_for_points(vector[:, :3], scale_parameter)
        corners = np.einsum('ikj,ij->ik', tetrahedron - vector[:, None, :3], vector[:, 3:])

        tetrahedron = self.apply_to_points(tetrahedron.reshape(-1, 3)).reshape(-1, 4, 3)
        v = gradient_from_tetrahedron(tetrahedron, corners)

        return v

    def add_abutting_fault(self, abutting_fault_feature, positive=None):

        # check whether the fault is on the hanging wall or footwall of abutting fault
        abutting_region = None
        if positive is None:
            pts = (
                self.__getitem__(0).builder.data[["X", "Y", "Z"]].to_numpy()
            )  # get_value_constraints()
            abut_value = np.nanmedian(abutting_fault_feature.evaluate_value(pts))
            positive = abut_value > 0
        # we want to crop the fault by the abutting
        # fault so create a positive/neg region and
        # include the fault centre and normal vector to help
        # outside of the fault interpolation support
        if positive:
            abutting_region = PositiveRegion(
                abutting_fault_feature,
                vector=abutting_fault_feature.fault_normal_vector,
                point=abutting_fault_feature.fault_centre,
            )
        if not positive:
            abutting_region = NegativeRegion(
                abutting_fault_feature,
                vector=abutting_fault_feature.fault_normal_vector,
                point=abutting_fault_feature.fault_centre,
            )
        self.abut[abutting_fault_feature.name] = abutting_region
        self.__getitem__(0).add_region(abutting_region)
        return abutting_region
    def save(self, filename, scalar_field=True, slip_vector=True, surface=True):
        """
        Save the fault to a file

        Parameters
        ----------
        filename - str
            filename to save to

        Returns
        -------

        """
        filename = str(filename)
        ext = filename.split(".")[-1]
        info = ''
        if scalar_field:
            if ext == '.geoh5':
                info = ''
            else:
                info = '_scalar_field'
            a = self.features[0].scalar_field()
            a.merge(self.features[1].scalar_field())
            a.merge(self.features[2].scalar_field())
            a.merge(self.displacementfeature.scalar_field())
            a.save(f'{filename}{info}.{ext}')
        if slip_vector:
            if ext == '.geoh5':
                info = ''
            else:
                info = '_slip_vector'
            self.vector_field().save(f'{filename}{info}.{ext}')
        if surface:
            if ext == '.geoh5':
                info = ''
            else:
                info = '_surface'
            self.surfaces([0])[0].save(f'{filename}{info}.{ext}')
