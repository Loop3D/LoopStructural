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
    def __init__(self, feature):
        self.feature = feature

    def __call__(self, xyz):
        return np.logical_or(
            self.feature.evaluate_value(xyz) > 0,
            np.isnan(self.feature.evaluate_value(xyz)),
        )


class NegativeRegion:
    def __init__(self, feature):
        self.feature = feature

    def __call__(self, xyz):
        return np.logical_or(
            self.feature.evaluate_value(xyz) < 0,
            np.isnan(self.feature.evaluate_value(xyz)),
        )
