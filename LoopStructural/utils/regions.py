import numpy as np


class BaseRegion:
    def __init__(self):
        self.type = "BaseRegion"


class RegionEverywhere(BaseRegion):
    def __init__(self):
        super().__init__()
        self.type = "RegionEverywhere"

    def __call__(self, xyz):
        return np.ones(xyz.shape[0], dtype="bool")


class RegionFunction(BaseRegion):
    def __init__(self, function):
        super().__init__()
        self.function = function

    def __call__(self, xyz):
        return self.function(xyz)


class PositiveRegion:
    """Helper class for evaluating whether you are in the positive  region of a scalar field.
    If its outside of the support it will interpolate the average gradient at a point on the 0 isovalue
    and calculate the distance from this. Alternatively, a point and vector can be used to save computational time
    """

    def __init__(self, feature, vector=None, point=None):
        self.feature = feature
        self.vector = vector
        self.point = point

    def __call__(self, xyz):
        val = self.feature.evaluate_value(xyz)
        # find a point on/near 0 isosurface
        if self.point is None:
            mask = np.zeros(xyz.shape[0], dtype="bool")
            mask[:] = val < 0
            if np.sum(mask) == 0:
                raise ValueError("Cannot find point on surface")
            centre = xyz[mask, :][0, :]
        else:
            centre = self.point
        if self.vector is None:
            average_gradient = self.feature.evaluate_gradient(np.array([centre]))[0]
            average_gradient[2] = 0
            average_gradient /= np.linalg.norm(average_gradient)
        else:
            average_gradient = self.vector
        # distance = ((xyz[:,0] - centre[None,0])*average_gradient[0] +
        #             (xyz[:,1] - centre[None,1])*average_gradient[1] +
        #             ( xyz[:,2] - centre[None,2])*average_gradient[2])
        distance = np.einsum("ij,j->i", centre[None, :] - xyz, average_gradient)
        return np.logical_or(
            np.logical_and(~np.isnan(val), val > 0),
            np.logical_and(np.isnan(val), distance > 0),
        )


class NegativeRegion:
    """Helper class for evaluating whether you are in the positive  region of a scalar field.
    If its outside of the support it will interpolate the average gradient at a point on the 0 isovalue
    and calculate the distance from this. Alternatively, a point and vector can be used to save computational time
    """

    def __init__(self, feature, vector=None, point=None):
        self.feature = feature
        self.vector = vector
        self.point = point

    def __call__(self, xyz):
        val = self.feature.evaluate_value(xyz)
        # find a point on/near 0 isosurface
        if self.point is None:
            mask = np.zeros(xyz.shape[0], dtype="bool")
            mask[:] = val < 0
            if np.sum(mask) == 0:
                raise ValueError("Cannot find point on surface")
            centre = xyz[mask, :][0, :]
        else:
            centre = self.point
        if self.vector is None:
            average_gradient = self.feature.evaluate_gradient(np.array([centre]))[0]
            average_gradient[2] = 0
            average_gradient /= np.linalg.norm(average_gradient)
        else:
            average_gradient = self.vector
        distance = np.einsum("ij,j->i", xyz - centre[None, :], average_gradient)
        # distance = ((xyz[:,0] - centre[None,0])*average_gradient[0] +
        #             (xyz[:,1] - centre[None,1])*average_gradient[1] +
        #             ( xyz[:,2] - centre[None,2])*average_gradient[2])
        # return np.logical_or(np.logical_and(~np.isnan(val),val
        #      < 0),
        #      np.logical_and(np.isnan(val),distance>0))
        return np.logical_or(
            np.logical_and(~np.isnan(val), val < 0),
            np.logical_and(np.isnan(val), distance < 0),
        )
