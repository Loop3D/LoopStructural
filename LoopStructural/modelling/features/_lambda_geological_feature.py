"""
Geological features
"""
from LoopStructural.modelling.features import BaseFeature
from LoopStructural.utils import getLogger
from LoopStructural.modelling.features import FeatureType
import numpy as np

logger = getLogger(__name__)


class LambdaGeologicalFeature(BaseFeature):
    def __init__(
        self,
        function=None,
        name="unnamed_lambda",
        gradient_function=None,
        model=None,
        regions=[],
        faults=[],
        builder=None,
    ):
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.type = FeatureType.LAMBDA
        self.function = function
        self.gradient_function = gradient_function

    def evaluate_value(self, xyz):
        v = np.zeros((xyz.shape[0]))
        if self.function is None:
            v[:] = np.nan
        else:
            v[:] = self.function(xyz)
        return v

    def evaluate_gradient(self, xyz):
        v = np.zeros((xyz.shape[0], 3))
        if self.gradient_function is None:
            v[:, :] = np.nan
        else:
            v[:, :] = self.gradient_function(xyz)
        return v
