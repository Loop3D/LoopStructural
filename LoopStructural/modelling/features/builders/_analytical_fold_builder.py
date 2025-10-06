from ._base_builder import BaseBuilder
from .._lambda_geological_feature import LambdaGeologicalFeature
import numpy as np
class AnalyticalFoldBuilder(BaseBuilder):
    def __init__(self, model, name: str = 'Feature'):
        super().__init__(model=model,name=name)
        self._wavelength = np.max(model.bounding_box.length)
        self._amplitude = np.min(model.bounding_box.length)
        self._centre = model.bounding_box
    @property
    def amplitude(self):
        return self._amplitude
    @property
    def wavelength(self):
        return self._wavelength
    
    @property
    def feature(self):
        def function(xyz):
            return xyz[:,2]+np.sin(xyz[:,0]/self.wavelength)*self.amplitude
        return LambdaGeologicalFeature(function=function,model=self.model,name=self.name)
