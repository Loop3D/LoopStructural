import numpy as np

class RegionEverywhere:
    def __call__(self,xyz):
        return np.ones(xyz.shape[0],dtype='bool')

class RegionFunction:
    def __init__(self,function):
        self.function = function
    def __call__(self,xyz):
        return self.function(xyz)