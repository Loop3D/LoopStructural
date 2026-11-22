from abc import ABC, abstractmethod

import numpy as np

from ._fault_block import block_kwargs, block_side, force_side


class BaseRegion(ABC):
    @abstractmethod
    def __init__(self, feature, vector=None, point=None):
        self.feature = feature
        self.vector = vector
        self.point = point
        self.name = None
        self.parent = None

    @abstractmethod
    def __call__(self, xyz, *, block=None) -> np.ndarray:
        """Evaluate the region based on the input coordinates.

        ``block`` is a fault block (see ``LoopStructural.utils._fault_block``);
        a region that depends on a faulted feature evaluates it in that block.
        """


class RegionEverywhere(BaseRegion):
    def __init__(self):
        super().__init__()
        self.type = "RegionEverywhere"

    def __call__(self, xyz, *, block=None):
        return np.ones(xyz.shape[0], dtype="bool")


class RegionFunction(BaseRegion):
    def __init__(self, function):
        super().__init__()
        self.function = function

    def __call__(self, xyz, *, block=None):
        return self.function(xyz)

class BaseSignRegion(BaseRegion):
    """Helper class for evaluating whether you are in the positive  region of a scalar field.
    If its outside of the support it will interpolate the average gradient at a point on the 0 isovalue
    and calculate the distance from this. Alternatively, a point and vector can be used to save computational time
    """

    def __init__(self, feature, vector=None, point=None):
        super().__init__(feature, vector, point)
        self.name = 'PositiveRegion'
        self.parent = feature

    def _calculate_value_and_distance(
        self, xyz, precomputed_val=None, block=None
    ) -> tuple[np.ndarray, np.ndarray]:
        if precomputed_val is not None:
            val = precomputed_val
        else:
            val = self.feature.evaluate_value(xyz, **block_kwargs(block))
        # find a point on/near 0 isosurface — compute once and cache on self
        if self.point is None:
            mask = np.zeros(xyz.shape[0], dtype="bool")
            mask[:] = val < 0
            if np.sum(mask) == 0:
                raise ValueError("Cannot find point on surface")
            self.point = xyz[mask, :][0, :]  # cache: characterises the fault, not the query points
        centre = self.point
        if self.vector is None:
            average_gradient = self.feature.evaluate_gradient(np.array([centre]))[0]
            average_gradient[2] = 0
            average_gradient /= np.linalg.norm(average_gradient)
            self.vector = average_gradient  # cache: fault geometry doesn't change between calls
        average_gradient = self.vector

        distance = np.einsum(
            "ij,j->i", centre[None, :] - xyz, average_gradient.reshape(-1, 3)[0, :]
        )
        # Where the feature has no support the distance decides the side, so a
        # fault block that fixes the side of this fault fixes it here too.
        distance = force_side(distance, block_side(block, getattr(self.feature, "name", None)))
        return val, distance


class PositiveRegion(BaseSignRegion):
    """Helper class for evaluating whether you are in the positive  region of a scalar field.
    If its outside of the support it will interpolate the average gradient at a point on the 0 isovalue
    and calculate the distance from this. Alternatively, a point and vector can be used to save computational time
    """

    def __init__(self, feature, vector=None, point=None):
        super().__init__(feature, vector, point)
        self.name = 'PositiveRegion'
        self.parent = feature

    def __call__(self, xyz, precomputed_val=None, *, block=None) -> np.ndarray:
        val, distance = self._calculate_value_and_distance(
            xyz, precomputed_val=precomputed_val, block=block
        )
        return np.logical_or(
            np.logical_and(~np.isnan(val), val > 0),
            np.logical_and(np.isnan(val), distance > 0),
        )


class NegativeRegion(BaseSignRegion):
    """Helper class for evaluating whether you are in the positive  region of a scalar field.
    If its outside of the support it will interpolate the average gradient at a point on the 0 isovalue
    and calculate the distance from this. Alternatively, a point and vector can be used to save computational time
    """

    def __init__(self, feature, vector=None, point=None):
        super().__init__(feature, vector, point)
        self.name = 'NegativeRegion'
        self.parent = feature

    def __call__(self, xyz, precomputed_val=None, *, block=None) -> np.ndarray:
        val, distance = self._calculate_value_and_distance(
            xyz, precomputed_val=precomputed_val, block=block
        )
        return np.logical_or(
            np.logical_and(~np.isnan(val), val < 0),
            np.logical_and(np.isnan(val), distance < 0),
        )
