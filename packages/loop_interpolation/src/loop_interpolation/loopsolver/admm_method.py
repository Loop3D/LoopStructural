import numpy as np


class ADMM:
    def __init__(self, nelements: int):
        if nelements < 1:
            raise ValueError("nelements must be greater than 0")
        self.z = np.zeros(nelements)
        self.u = np.zeros(nelements)
        self.nelements = nelements

    def admm_method_iterate_admm_array(self, xmin, xmax, x):

        arg = np.nan
        inside = False
        arg = x[:, 0] + self.u
        inside = np.logical_and(arg >= xmin[:, 0], arg <= xmax[:, 0])
        self.z[inside] = arg[inside]
        below = np.logical_and(arg < xmin[:, 0], ~inside)
        above = np.logical_and(arg > xmax[:, 0], ~inside)
        self.z[below] = xmin[below, 0]
        self.z[above] = xmax[above, 0]
        # calculate the
        self.u = self.u + x[:, 0] - self.z
        return self.z - self.u
