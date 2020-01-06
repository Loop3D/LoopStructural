import logging

from LoopStructural.modelling.fault.fault_function_feature import FaultDisplacementFeature

logger = logging.getLogger(__name__)
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
        A slip event of a faut
        Parameters
        ----------
        faultframe - the fault frame defining the faut geometry
        faultfunction :
        steps : int
        kwargs
        """
        self.faultframe = faultframe
        self.type = 'fault'
        self.name = kwargs.get('name', self.faultframe.name)
        self.displacement = displacement
        self.faultfunction = faultfunction
        self.steps = steps
        self.regions = []

        self.displacementfeature = None
        if self.faultframe is not None:
            self.displacementfeature = FaultDisplacementFeature(
                self.faultframe, self.faultfunction, name = self.name)

    def __getitem__(self, item):
        return self.faultframe[item]

    def add_region(self, region):
        """

        Parameters
        ----------
        region - boolean function(x,y,z)
                - returns true if inside region, false if outside
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
        return self.faultframe[0].evaluate_value(locations[mask, :])

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
        v = np.zeros(locations.shape[0])
        v[:] = np.nan
        mask = np.zeros(locations.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask, r(locations))
        # need to scale with fault displacement
        return self.faultframe[1].evaluate_gradient(locations[mask, :])

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
        d[np.isnan(gx)] = 0
        # d[~np.isnan(gx)][gx[~np.isnan(gx)]>0] = 1
        d[gx > 0] = 1.
        d[gx < 0] = 0.
        if self.faultfunction is not None:
            d = self.faultfunction(gx, gy, gz)
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
            d[np.isnan(gx)] = 0
            d[gx > 0] = 1.
            d[gx < 0] = 0.
            if self.faultfunction is not None:
                d = self.faultfunction(gx, gy, gz)
            d *= self.displacement
            # normalise when length is >0
            g_mag = np.linalg.norm(g, axis=1)
            g[g_mag > 0.] /= g_mag[g_mag > 0, None]
            # multiply displacement vector by the displacement magnitude for
            # step
            g *= (1. / steps) * d[:, None]

            # apply displacement
            newp[mask, :] += g
        return newp

    def apply_to_data(self, data):
        """
        Unfault the data in the list provided

        Parameters
        ----------
        data - list containing Points

        Returns
        -------

        """
        logger.info("Applying fault")
        steps = self.steps
        # TODO make this numpy arrays
        if data is None:
            return
        for d in data:
            hw = self.faultframe.features[0].evaluate_value(
                np.array([d.pos])) > 0
            for i in range(steps):
                g = self.faultframe.features[1].evaluate_gradient(
                    np.array([d.pos]))
                length = np.linalg.norm(g, axis=1)
                g[np.logical_and(~np.isnan(length), length > 0), :] /= length[
                    np.logical_and(~np.isnan(length), length > 0), None]
                if self.faultfunction is None and hw:
                    g *= (1. / steps) * 1. * self.displacement
                    g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
                    d.pos = d.pos + g[0]
