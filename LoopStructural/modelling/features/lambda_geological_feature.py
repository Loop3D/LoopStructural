"""
Geological features
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class LambdaGeologicalFeature:
    def __init__(self,function = None,name = 'unnamed_lambda', gradient_function = None, model = None):
        self.function = function
        self.name = name
        self.gradient_function = gradient_function
        self.model = model
        if self.model is not None:
            v = function(self.model.regular_grid((10, 10, 10)))
            self._min = np.nanmin(v)#function(self.model.regular_grid((10, 10, 10))))
            self._max = np.nanmax(v)
        else:
            self._min = 0
            self._max = 0
    def evaluate_value(self, xyz):
        v = np.zeros((xyz.shape[0]))
        if self.function is None:
            v[:] =  np.nan
        else:
            v[:] = self.function(xyz)
        return v
    def evaluate_gradient(self,xyz):
        v = np.zeros((xyz.shape[0],3))
        if self.gradient_function is None:
            v[:,:] =  np.nan
        else:
            v[:,:] = self.gradient_function(xyz)
        return v
    def min(self):
        return self._min
    def max(self):
        return self._max

    