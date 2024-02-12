from ....modelling.features.fault._fault_function_feature import (
    FaultDisplacementFeature,
)
from ....modelling.features import FeatureType
from ....modelling.features.fault._fault_function import BaseFault, BaseFault3D
from ....utils import getLogger, NegativeRegion, PositiveRegion
from ....modelling.features import StructuralFrame

from concurrent.futures import ThreadPoolExecutor
import numpy as np

logger = getLogger(__name__)

use_threads = True


class FaultSegment(StructuralFrame):
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, features, name, faultfunction=None, steps=10, displacement=1.0, fold=None):
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
        StructuralFrame.__init__(self, features, name, fold)
        self.type = FeatureType.FAULT
        self.displacement = displacement
        self._faultfunction = BaseFault.fault_displacement
        self.steps = steps
        self.regions = []
        self.faults_enabled = True

        self.builder = None
        self.splay = {}
        self.abut = {}

    @property
    def faultfunction(self):
        return self._faultfunction

    @faultfunction.setter
    def faultfunction(self, value):
        if callable(value):
            self._faultfunction = value
        elif isinstance(value, str) and value == "BaseFault":
            self._faultfunction = BaseFault.fault_displacement
        elif isinstance(value, str) and value == "BaseFault3D":
            self._faultfunction = BaseFault3D.fault_displacement
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
        return FaultDisplacementFeature(self, self.faultfunction, name=self.name)

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

    def evaluate_value(self, locations):
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
        return self.__getitem__(0).evaluate_value(locations)

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
            except:
                logger.error("nan slicing ")
        # need to scale with fault displacement
        v[mask, :] = self.__getitem__(0).evaluate_gradient(locations[mask, :])
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
            d[mask] = self.faultfunction(gx[mask], gy[mask], gz[mask])
        return d * self.displacement

    def apply_to_points(self, points, reverse=False):
        """
        Unfault the array of points

        Parameters
        ----------
        points - numpy array Nx3

        Returns
        -------

        """
        steps = self.steps
        newp = np.copy(points).astype(float)
        # evaluate fault function for all points
        # then define mask for only points affected by fault
        gx = None
        gy = None
        gz = None
        if use_threads:
            with ThreadPoolExecutor(max_workers=8) as executor:
                # all of these operations should be
                # independent so just run as different threads
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
            d[mask] = self.faultfunction(gx[mask], gy[mask], gz[mask])
        mask = np.abs(d) > 0.0

        d *= self.displacement
        if reverse:
            d *= -1.0
        # calculate the fault frame for the evaluation points
        for _i in range(steps):
            gx = None
            gy = None
            gz = None
            g = None
            if use_threads:
                with ThreadPoolExecutor(max_workers=8) as executor:
                    # all of these operations should be
                    # independent so just run as different threads
                    gx_future = executor.submit(self.__getitem__(0).evaluate_value, newp[mask, :])
                    g_future = executor.submit(self.__getitem__(1).evaluate_gradient, newp[mask, :])
                    gy_future = executor.submit(self.__getitem__(1).evaluate_value, newp[mask, :])
                    gz_future = executor.submit(self.__getitem__(2).evaluate_value, newp[mask, :])
                    gx = gx_future.result()
                    g = g_future.result()
                    gy = gy_future.result()
                    gz = gz_future.result()
            else:
                gx = self.__getitem__(0).evaluate_value(newp[mask, :])
                gy = self.__getitem__(1).evaluate_value(newp[mask, :])
                gz = self.__getitem__(2).evaluate_value(newp[mask, :])
                g = self.__getitem__(1).evaluate_gradient(newp[mask, :])
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

        regular_tetrahedron = np.array(
            [
                [np.sqrt(8 / 9), 0, -1 / 3],
                [-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3],
                [-np.sqrt(2 / 9), -np.sqrt(2 / 3), -1 / 3],
                [0, 0, 1],
            ]
        )
        regular_tetrahedron *= scale_parameter
        xyz = vector[:, :3]
        tetrahedron = np.zeros((xyz.shape[0], 4, 3))
        tetrahedron[:] = xyz[:, None, :]
        tetrahedron[:, :, :] += regular_tetrahedron[None, :, :]

        vectors = vector[:, 3:]
        corners = np.einsum('ikj,ij->ik', tetrahedron - xyz[:, None, :], vectors)
        tetrahedron = tetrahedron.reshape(-1, 3)
        tetrahedron = self.apply_to_points(tetrahedron)
        tetrahedron = tetrahedron.reshape(-1, 4, 3)
        m = np.array(
            [
                [
                    (tetrahedron[:, 1, 0] - tetrahedron[:, 0, 0]),
                    (tetrahedron[:, 1, 1] - tetrahedron[:, 0, 1]),
                    (tetrahedron[:, 1, 2] - tetrahedron[:, 0, 2]),
                ],
                [
                    (tetrahedron[:, 2, 0] - tetrahedron[:, 0, 0]),
                    (tetrahedron[:, 2, 1] - tetrahedron[:, 0, 1]),
                    (tetrahedron[:, 2, 2] - tetrahedron[:, 0, 2]),
                ],
                [
                    (tetrahedron[:, 3, 0] - tetrahedron[:, 0, 0]),
                    (tetrahedron[:, 3, 1] - tetrahedron[:, 0, 1]),
                    (tetrahedron[:, 3, 2] - tetrahedron[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        v = np.sum(element_gradients * corners, axis=2)
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
