from concurrent.futures import ThreadPoolExecutor

import numpy as np

from LoopStructural.utils.maths import gradient_from_tetrahedron, regular_tetraherdron_for_points

from ....modelling.features import FeatureType, StructuralFrame
from ....modelling.features.fault._fault_function import (
    BaseFault,
    BaseFault3D,
    FaultDisplacement,
    block_displacement,
)
from ....modelling.features.fault._fault_function_feature import (
    FaultDisplacementFeature,
)
from ....utils import NegativeRegion, PositiveRegion, getLogger, public_api
from ....utils._fault_block import block_kwargs, block_side, check_block, force_side

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
    def faultfunction(self):
        return self._faultfunction

    @faultfunction.setter
    def faultfunction(self, value):
        if isinstance(value, str) and value == "BaseFault":
            self._faultfunction = BaseFault().fault_displacement
        elif isinstance(value, str) and value == "BaseFault3D":
            self._faultfunction = BaseFault3D().fault_displacement
        elif callable(value):
            self._faultfunction = value
        else:
            raise ValueError("Fault function must be a function or BaseFault")

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

            # The fault centre is a world point. Do not convert it again.
            pts = self.fault_centre[None, :]
            fault_ellipsoid = pv.PolyData(pts)
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
        elif not scale:
            self.displacement = displacement
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

    def evaluate_value(self, locations, ignore_regions=False, *, block=None):
        """
        Return the value of the fault surface scalar field

        Parameters
        ----------
        locations - numpy array
            location to evaluate scalar field
        ignore_regions : bool
            do not apply the regions of the fault (for example an abutting fault)
        block : dict, optional
            fault block, maps fault names to +1 (hanging wall), -1 (footwall) or
            0 (absent); see :meth:`GeologicalModel.evaluate_feature_block_value`.
            If the block names this fault with +1 or -1, the value is put on that
            side (``side * |gx|``), so a region defined by this fault gives the
            same answer at every point.

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
        block = check_block(block)
        v = self.__getitem__(0).evaluate_value(
            locations, ignore_regions=ignore_regions, **block_kwargs(block)
        )
        return force_side(v, block_side(block, self.name))

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

    def _evaluate_frame(self, points, block=None, precomputed_gx=None, restored_through=()):
        """Evaluate the three fault frame coordinates at the points.

        Parameters
        ----------
        points : np.ndarray
            Nx3 locations
        block : dict, optional
            checked fault block; gx is put on the side that the block gives
            for this fault, and the frame is evaluated in the block of the
            faults that displace this fault
        precomputed_gx : np.ndarray, optional
            gx already evaluated at the same points
        restored_through : sequence of FaultSegment
            faults whose restoration the points already have; the frame does
            not apply them again

        Returns
        -------
        tuple of np.ndarray
            gx, gy, gz
        """
        kwargs = block_kwargs(block)
        if restored_through:
            kwargs["restored_through"] = tuple(restored_through)
        if precomputed_gx is not None:
            gx = precomputed_gx
            gy = self.__getitem__(1).evaluate_value(points, **kwargs)
            gz = self.__getitem__(2).evaluate_value(points, **kwargs)
        elif use_threads:
            with ThreadPoolExecutor(max_workers=8) as executor:
                # all of these operations should be independent so
                # just run as different threads
                gx_future = executor.submit(self.__getitem__(0).evaluate_value, points, **kwargs)
                gy_future = executor.submit(self.__getitem__(1).evaluate_value, points, **kwargs)
                gz_future = executor.submit(self.__getitem__(2).evaluate_value, points, **kwargs)
                gx = gx_future.result()
                gy = gy_future.result()
                gz = gz_future.result()
        else:
            gx = self.__getitem__(0).evaluate_value(points, **kwargs)
            gy = self.__getitem__(1).evaluate_value(points, **kwargs)
            gz = self.__getitem__(2).evaluate_value(points, **kwargs)
        return gx, gy, gz

    def _displacement_fraction(self, gx, gy, gz, block=None, offset=0.0):
        """The displacement at the frame coordinates, as a fraction of the displacement.

        Parameters
        ----------
        gx, gy, gz : np.ndarray
            fault frame coordinates (:meth:`_evaluate_frame`), NaN where the
            fault has no frame
        block : dict, optional
            checked fault block. If it gives a side for this fault, every
            point takes the displacement of that side: the profile itself on
            that side and its smooth continuation across the fault on the
            other side (see ``_fault_function.block_displacement``).
        offset : float
            added to gx in the profile

        Returns
        -------
        np.ndarray
            zero where a frame coordinate is NaN
        """
        side = block_side(block, self.name)
        forced = force_side(gx, side) if side in (-1, 1) else gx
        d = np.zeros(gx.shape)
        mask = np.logical_and(~np.isnan(forced), ~np.isnan(gy))
        mask = np.logical_and(mask, ~np.isnan(gz))
        d[~mask] = 0
        gx_mask = np.zeros_like(mask, dtype=bool)
        gx_mask[mask] = forced[mask] > 0
        d[gx_mask] = 1.0
        if self.faultfunction is not None:
            d[mask] = block_displacement(
                self.faultfunction, gx[mask] + offset, forced[mask] + offset, gy[mask], gz[mask],
                side,
            )
        return d

    def evaluate_displacement(self, points, *, block=None):
        """Evaluate the displacement magnitude of the fault at the points.

        Parameters
        ----------
        points : np.ndarray
            Nx3 locations
        block : dict, optional
            fault block (see :meth:`apply_to_points`)

        Returns
        -------
        np.ndarray
            displacement magnitude, zero where the fault does not move points
        """
        block = check_block(block)
        if block_side(block, self.name) == 0:
            return np.zeros(np.shape(points)[0])
        newp = np.copy(points).astype(float)
        # evaluate fault function for all
        # points then define mask for only points affected by fault
        gx, gy, gz = self._evaluate_frame(newp, block)
        d = self._displacement_fraction(gx, gy, gz, block, offset=self.fault_offset)
        return d * self.displacement

    @public_api(tier="provisional")
    def evaluate_frame_value(self, points, coordinate: int = 0, *, block=None) -> np.ndarray:
        """Evaluate one coordinate of the fault frame, without the fault regions.

        Coordinate 0 (gx) is zero on the fault surface and positive in the
        hanging wall, so ``gx / |grad gx|`` is a signed distance to the fault
        near the fault. Coordinate 1 (gy) follows the slip direction and
        coordinate 2 (gz) the fault extent. The regions of the fault (for
        example the side of an abutting fault where the fault exists) are not
        applied: the fault surface is extended over the whole support.

        Parameters
        ----------
        points : np.ndarray
            Nx3 locations, world coordinates
        coordinate : int
            0, 1 or 2
        block : dict, optional
            fault block of the younger faults that displace this fault (see
            :meth:`GeologicalModel.evaluate_feature_block_value`). The frame
            is then the smooth frame of that block, with no jump across the
            younger faults. A side given for this fault itself has no effect.

        Returns
        -------
        np.ndarray
            N values, NaN where the frame has no support
        """
        block = check_block(block)
        points = np.asarray(points, dtype=float)
        return self.__getitem__(coordinate).evaluate_value(
            points, ignore_regions=True, **block_kwargs(block)
        )

    @public_api(tier="provisional")
    def evaluate_frame_gradient(self, points, coordinate: int = 0, *, block=None) -> np.ndarray:
        """Evaluate the gradient of one coordinate of the fault frame.

        See :meth:`evaluate_frame_value` for the parameters.

        Returns
        -------
        np.ndarray
            Nx3 gradients, NaN where the frame has no support
        """
        block = check_block(block)
        points = np.asarray(points, dtype=float)
        return self.__getitem__(coordinate).evaluate_gradient(
            points, ignore_regions=True, **block_kwargs(block)
        )

    def _tip_ends(self) -> dict:
        """Find the ends of gy and gz where the displacement profile is zero.

        Returns
        -------
        dict of int to list of int
            For coordinate 1 and 2, the ends (+1 and/or -1) where the profile
            is zero at |g| >= 1. An empty list means no tip on that coordinate.
        """
        if isinstance(self.faultfunction, FaultDisplacement):
            return {c: list(e) for c, e in self.faultfunction.tip_ends.items()}
        # A plain function of (gx, gy, gz): its tip is not known.
        return {1: [], 2: []}

    def _evaluate_tip(self, points, block, gradient: bool):
        """Evaluate the tip function and, if asked, its gradient."""
        points = np.asarray(points, dtype=float)
        value = np.full(points.shape[0], np.inf)
        grad = np.zeros(points.shape)
        for coordinate, ends in self._tip_ends().items():
            if not ends:
                continue
            g = self.evaluate_frame_value(points, coordinate, block=block)
            if len(ends) == 2:
                part, sign = 1.0 - np.abs(g), -np.sign(g)
            else:
                part, sign = 1.0 - ends[0] * g, -ends[0] * np.ones_like(g)
            if gradient:
                dg = self.evaluate_frame_gradient(points, coordinate, block=block)
                lower = (part < value) | np.isnan(part)
                grad[lower] = sign[lower, None] * dg[lower]
            value = np.where(np.isnan(part) | (part < value), part, value)
        return value, grad

    @public_api(tier="provisional")
    def evaluate_tip_value(self, points, *, block=None) -> np.ndarray:
        """Evaluate the tip function of a finite fault.

        The function is ``min(1 - |gy|, 1 - |gz|)`` over the frame coordinates
        whose displacement profile goes to zero at ``|g| = 1`` (``1 - g`` or
        ``1 + g`` if the profile is zero at one end only). It is positive where
        the fault displaces points, zero on the tip line and negative beyond.
        Beyond the tip the displacement is exactly zero, so the fields of the
        two fault blocks are equal there.

        Parameters
        ----------
        points : np.ndarray
            Nx3 locations, world coordinates
        block : dict, optional
            fault block of the faults that displace this fault (see
            :meth:`evaluate_frame_value`)

        Returns
        -------
        np.ndarray
            N values; ``inf`` everywhere if the fault has no tip (its profile
            is not zero at ``|g| = 1``, or the fault function is not a
            :class:`FaultDisplacement`); NaN where the frame has no support
        """
        return self._evaluate_tip(points, check_block(block), gradient=False)[0]

    @public_api(tier="provisional")
    def evaluate_tip_gradient(self, points, *, block=None) -> np.ndarray:
        """Evaluate the gradient of :meth:`evaluate_tip_value`.

        Returns
        -------
        np.ndarray
            Nx3 gradients of the frame coordinate that gives the minimum;
            zero if the fault has no tip
        """
        return self._evaluate_tip(points, check_block(block), gradient=True)[1]

    def apply_to_points(
        self, points, reverse=False, precomputed_gx=None, *, block=None, restored_through=()
    ):
        """
        Unfault the array of points

        Parameters
        ----------
        points - numpy array Nx3
        precomputed_gx - optional pre-evaluated gx values (same points, avoids duplicate evaluation)
        block : dict, optional
            fault block: maps fault names to +1 (hanging wall), -1 (footwall)
            or 0 (absent). If it names this fault, every point is restored as
            if it were on that side (0: the points are not moved). A point on
            the other side takes the smooth continuation of the displacement
            profile of that side across the fault, so the restored field is
            C2 across the fault surface (see
            ``_fault_function.block_displacement``). The block
            is also used to evaluate the frame of this fault, so a fault that
            younger faults displace is evaluated in the block of those faults.
            The call changes no state.
        restored_through : sequence of FaultSegment, optional
            faults whose restoration the points already have. The frame of
            this fault was built in the space restored by the younger faults
            that displace it, so it is evaluated at the points as they are,
            without applying those faults a second time.

        Returns
        -------

        """
        logger.info(f'Applying fault {self.name} to points {points.shape}')
        block = check_block(block)
        if block_side(block, self.name) == 0:
            return np.copy(points).astype(float)
        kwargs = block_kwargs(block)
        if restored_through:
            kwargs["restored_through"] = tuple(restored_through)
        steps = self.steps
        newp = np.copy(points).astype(float)
        # evaluate fault function for all points
        # then define mask for only points affected by fault
        # gx may be supplied by caller to avoid re-evaluation (precomputed from region check)
        gx, gy, gz = self._evaluate_frame(
            newp, block, precomputed_gx=precomputed_gx, restored_through=restored_through
        )
        d = self._displacement_fraction(gx, gy, gz, block)
        mask = np.abs(d) > 0.0

        d *= self.displacement
        if reverse:
            d *= -1.0
        # calculate the fault frame for the evaluation points
        for _i in range(steps):
            g = None
            if _i == 0:
                # Reuse gx/gy/gz from the initial full-array evaluation above — newp[mask] hasn't
                # moved yet on the first iteration, so values are identical.
                gx_m = gx[mask]
                gy_m = gy[mask]
                gz_m = gz[mask]
            else:
                gx_m, gy_m, gz_m = self._evaluate_frame(
                    newp[mask, :], block, restored_through=restored_through
                )
            g = self.__getitem__(1).evaluate_gradient(newp[mask, :], ignore_regions=True, **kwargs)
            gx, gy, gz = gx_m, gy_m, gz_m  # alias for block below
            # # get the fault frame val/grad for the points
            # determine displacement magnitude, for constant displacement
            # hanging wall should be > 0
            mask2 = np.logical_and(~np.isnan(gx), ~np.isnan(gy))
            mask2 = np.logical_and(mask2, ~np.isnan(gz))
            # the gradient is NaN for points on a younger fault: do not move those points
            mask2 &= np.all(np.isfinite(g), axis=1)
            g[~mask2] = 0.0
            d = self._displacement_fraction(gx, gy, gz, block)
            d *= self.displacement
            if reverse:
                # forward (geological) movement: the opposite of the restoration
                d *= -1.0
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

    def apply_to_vectors(
        self, vector: np.ndarray, scale_parameter: float = 1.0, *, restored_through=()
    ) -> np.ndarray:
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

        kwargs = {"restored_through": tuple(restored_through)} if restored_through else {}
        tetrahedron = self.apply_to_points(tetrahedron.reshape(-1, 3), **kwargs).reshape(-1, 4, 3)
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
