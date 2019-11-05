import numpy as np
import logging
logger = logging.getLogger(__name__)


class FaultSegment:
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, faultframe, **kwargs):  # mesh,fault_event,data,name,region):
        """
        A slip event of a faut
        Parameters
        ----------
        faultframe - the fault frame defining the faut geometry
        kwargs -
        displacement magnitude of displacement
        faultfunction F(fault_frame) describes displacement inside fault frame
        steps - how many steps to apply the fault on
        """
        self.faultframe = faultframe
        self.type = 'fault'
        self.name = kwargs.get('name',self.faultframe.name)
        self.displacement = kwargs.get('displacement', 1.)
        self.gy_min = kwargs.get('gy_min', -9999)
        self.gy_max = kwargs.get('gy_max', 9999)
        self.faultfunction = kwargs.get('faultfunction', None)
        self.steps = kwargs.get('steps', 10)
        self.regions = []

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
        print('eval')
        for r in self.regions:
            print('regio')
            mask = np.logical_and(mask,r(locations))
        return self.faultframe[0].evaluate_value(locations[mask,:])

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
            mask = np.logical_and(mask,r(locations))
        self.faultframe[1].evaluate_gradient(locations[mask,:])

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
        # calculate the fault frame for the evaluation points
        for i in range(steps):
            # get the fault frame val/grad for the points
            gx = self.faultframe.features[0].evaluate_value(points)
            g = self.faultframe.features[1].evaluate_gradient(points)
            gy = self.faultframe.features[1].evaluate_value(points)
            gz = self.faultframe.features[2].evaluate_value(points)
            # determine displacement magnitude, for constant displacement
            # hanging wall should be > 0

            d = np.zeros(gx.shape)
            d[np.isnan(gx)] = 0
            # d[~np.isnan(gx)][gx[~np.isnan(gx)]>0] = 1
            d[gx > 0] = 1.
            if self.faultfunction is not None:
                d = self.faultfunction(gx, gy, gz)
            d *= self.displacement

            # normalise when length is >0
            g_mag = np.linalg.norm(g, axis=1)
            g[g_mag>0.] /= g_mag[g_mag >0,None]
            # multiply displacement vector by the displacement magnitude for step
            g *= (1./steps)*d[:,None]
            
            # apply displacement
            newp += g
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
        steps = self.steps
        # TODO make this numpy arrays
        if data is None:
            return
        for d in data:
            hw = self.faultframe.features[0].evaluate_value(np.array([d.pos])) > 0
            for i in range(steps):
                g = self.faultframe.features[1].evaluate_gradient(np.array([d.pos]))
                length = np.linalg.norm(g, axis=1)
                g[np.logical_and(~np.isnan(length),length>0),:] /= length[np.logical_and(~np.isnan(length),length>0),None]
                if self.faultfunction is None and hw:
                    g *= (1. / steps) * 1.*self.displacement
                    g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
                    d.pos = d.pos + g[0]

    def slice(self, isovalue, bounding_box = None, nsteps = None):
        return self.faultframe[0].slice(isovalue, bounding_box = bounding_box, nsteps = nsteps)