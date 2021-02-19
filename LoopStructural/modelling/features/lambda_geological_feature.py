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
        return 0
    def max(self):
        return 0

    