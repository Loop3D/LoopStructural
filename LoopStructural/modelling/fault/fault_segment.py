import logging

from LoopStructural.modelling.fault.fault_function_feature import FaultDisplacementFeature
from LoopStructural.modelling.fault.fault_function import BaseFault
from LoopStructural.utils import getLogger
logger = getLogger(__name__)
from concurrent.futures import ThreadPoolExecutor
import numpy as np


class FaultSegment:
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, faultframe,
                 faultfunction = None,
                 steps = 3,
                 displacement=1.,
                 **kwargs):
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
        self.faultframe = faultframe
        self.type = 'fault'
        self.name = kwargs.get('name', self.faultframe.name)
        self.displacement = displacement
        self.faultfunction = faultfunction
        if faultfunction == 'BaseFault':
            self.faultfunction = BaseFault.fault_displacement
        self.steps = steps
        self.regions = []
        self.faults_enabled = True
        self.displacementfeature = None
        self.model = None
        if self.faultframe is not None:
            self.displacementfeature = FaultDisplacementFeature(
                self.faultframe, self.faultfunction, name = self.name)
        self.builder = None
    
    def __getitem__(self, item):
        """

        Parameters
        ----------
        item

        Returns
        -------

        """
        return self.faultframe[item]
    
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

    def set_displacement(self, displacement, scale = True):
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

    def toggle_faults(self):
        """
        Toggle faults that affect this fault segment

        Returns
        -------

        """
        self.faults_enabled = ~self.faults_enabled
        for i in range(3):
            self.faultframe[i].toggle_faults()

    def add_region(self, region):
        """

        Parameters
        ----------
        region : boolean function(x,y,z)
                A function that returns true if inside region, false if outside
                can be passed as a lambda function e.g.
                lambda pos : feature.evaluate_value(pos) > 0
        Returns
        -------

        """

        self.regions.append(region)

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

        return self.faultframe.features[0].evaluate_value(locations) > 0

    def inside_volume(self,locations):
        v = self.faultframe.evaluate_value(locations)
        return np.all(np.logical_and(v > -1,v<1),axis=1)

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
        v = np.zeros(locations.shape[0])
        v[:] = np.nan
        mask = np.zeros(locations.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask, r(locations))
        v[mask]=self.faultframe[0].evaluate_value(locations[mask, :])
        return v

    def mean(self):
        return self.faultframe[0].mean()

    def max(self):
        return self.faultframe[0].max()

    def min(self):
        return self.faultframe[0].min()

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
        v[:,:] = np.nan
        mask = np.zeros(locations.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask, r(locations))
        # need to scale with fault displacement
        v[mask,:] = self.faultframe[1].evaluate_gradient(locations[mask, :])
        return v

    def evaluate_displacement(self, points):
        newp = np.copy(points).astype(float)
        # evaluate fault function for all points then define mask for only points affected by fault
        with ThreadPoolExecutor(max_workers=8) as executor:
            # all of these operations should be independent so just run as different threads
            gx_future = executor.submit(self.faultframe.features[0].evaluate_value, newp)
            gy_future = executor.submit(self.faultframe.features[1].evaluate_value, newp)
            gz_future = executor.submit(self.faultframe.features[2].evaluate_value, newp)
            gx = gx_future.result()
            gy = gy_future.result()
            gz = gz_future.result()
        d = np.zeros(gx.shape)
        mask = np.logical_and(~np.isnan(gx),~np.isnan(gy))
        mask = np.logical_and(mask,~np.isnan(gz))
        d[~mask] = 0
        gx_mask = np.zeros_like(mask,dtype=bool)
        gx_mask[mask] = gx[mask] > 0
        d[gx_mask] = 1.
        if self.faultfunction is not None:
            d[mask] = self.faultfunction(gx[mask], gy[mask], gz[mask])
        return d
    def apply_to_points(self, points):
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
        # evaluate fault function for all points then define mask for only points affected by fault
        with ThreadPoolExecutor(max_workers=8) as executor:
            # all of these operations should be independent so just run as different threads
            gx_future = executor.submit(self.faultframe.features[0].evaluate_value, newp)
            gy_future = executor.submit(self.faultframe.features[1].evaluate_value, newp)
            gz_future = executor.submit(self.faultframe.features[2].evaluate_value, newp)
            gx = gx_future.result()
            gy = gy_future.result()
            gz = gz_future.result()
        d = np.zeros(gx.shape)
        mask = np.logical_and(~np.isnan(gx),~np.isnan(gy))
        mask = np.logical_and(mask,~np.isnan(gz))
        d[~mask] = 0
        gx_mask = np.zeros_like(mask,dtype=bool)
        gx_mask[mask] = gx[mask] > 0
        d[gx_mask] = 1.
        if self.faultfunction is not None:
            d[mask] = self.faultfunction(gx[mask], gy[mask], gz[mask])
        mask = np.abs(d) > 0.

        d *= self.displacement
        # calculate the fault frame for the evaluation points
        for i in range(steps):
            with ThreadPoolExecutor(max_workers=8) as executor:
                # all of these operations should be independent so just run as different threads
                gx_future = executor.submit(self.faultframe.features[0].evaluate_value, newp[mask, :])
                g_future = executor.submit(self.faultframe.features[1].evaluate_gradient, newp[mask, :])
                gy_future = executor.submit(self.faultframe.features[1].evaluate_value, newp[mask, :])
                gz_future = executor.submit(self.faultframe.features[2].evaluate_value, newp[mask, :])
                gx = gx_future.result()
                g = g_future.result()
                gy = gy_future.result()
                gz = gz_future.result()
            # # get the fault frame val/grad for the points
            # determine displacement magnitude, for constant displacement
            # hanging wall should be > 0
            d = np.zeros(gx.shape)
            mask2 = np.logical_and(~np.isnan(gx), ~np.isnan(gy))
            mask2 = np.logical_and(mask2, ~np.isnan(gz))
            d[~mask2] = 0
            gx_mask2 = np.zeros_like(mask2,dtype=bool)
            gx_mask2[mask2] = gx[mask2] > 0
            # d[~np.isnan(gx)][gx[~np.isnan(gx)]>0] = 1
            d[gx_mask2] = 1.
            # d[mask2][gx[mask2] < 0] = 0.
            # d[gx < 0] = 0.
            if self.faultfunction is not None:
                d[mask2] = self.faultfunction(gx[mask2], gy[mask2], gz[mask2])
            d *= self.displacement
            # normalise when length is >0
            g_mag = np.zeros(g.shape[0])
            g_mag[mask2] = np.linalg.norm(g[mask2], axis=1)
            # g_mag = np.linalg.norm(g[mask2], axis=1)
            g[g_mag > 0.] /= g_mag[g_mag > 0, None]
            # multiply displacement vector by the displacement magnitude for
            # step
            g *= (1. / steps) * d[:, None]

            # apply displacement
            newp[mask, :] += g
        return newp
