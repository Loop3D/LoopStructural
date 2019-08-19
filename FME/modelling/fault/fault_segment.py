import numpy as np

class FaultSegment:
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, faultframe,**kwargs):  # mesh,fault_event,data,name,region):
        """
        constructor
        :param faultframe: StructuralFrame object with three coordinates
        :param kwargs:
        """
        self.faultframe = faultframe
        self.displacement = 1
        self.gy_max = 9999
        self.gy_min = -9999
        self.steps = 10
        self.faultfunction = None
        if 'name' in kwargs:
            self.name = kwargs['name']
        if 'faultfunction' in kwargs:
            self.faultfunction = kwargs['faultfunction']
        if 'displacement' in kwargs:
            self.displacement = kwargs['displacement']
        if 'gy_min' in kwargs:
            self.gy_min = kwargs['gy_min']
        if 'gy_max' in kwargs:
            self.gy_max = kwargs['gy_max']
        if 'steps' in kwargs:
            self.steps = kwargs['steps']

    def evaluate(self, locations):
        return self.faultframe.features[0].evaluate_value(locations) > 0

    def apply_to_points(self, points):
        steps = self.steps
        newp = np.copy(points).astype(float)
        # calculate the fault frame for the evaluation points
        for i in range(steps):
            gx = self.faultframe.features[0].evaluate_value(points)
            g = self.faultframe.features[1].evaluate_gradient(points)
            gy = self.faultframe.features[1].evaluate_value(points)
            gz = self.faultframe.features[2].evaluate_value(points)
            d = np.zeros(gx.shape)
            d[gx > 0] = 1.
            if self.faultfunction is not None:
                d = self.faultfunction(gx, gy, gz)
            d *= self.displacement
            # cast newpoints as a float because python types suck

            # calc length of displacement vector
            g_mag = np.linalg.norm(g, axis=1)
            # normalise when length is >0
            g[g_mag>0.] /= g_mag[g_mag >0,None]
            g *= (1./steps)*d[:,None]
            newp += g
        return newp

    def apply_to_data(self, data):
        steps = self.steps
        # TODO make this numpy arrays
        if data is None:
            return
        for d in data:
            hw = self.faultframe.features[0].evaluate_value(np.array([d.pos])) > 0
            for i in range(steps):
                g = self.faultframe.features[1].evaluate_gradient(np.array([d.pos]))
                g /= np.linalg.norm(g, axis=1)[:, None]
                if self.faultfunction is None and hw:
                    g *= (1. / steps) * 1.*self.displacement
                    g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
                    d.pos = d.pos + g[0]
        return
