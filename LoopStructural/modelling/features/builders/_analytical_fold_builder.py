import numpy as np

from .._lambda_geological_feature import LambdaGeologicalFeature
from ._base_builder import BaseBuilder


class AnalyticalFoldBuilder(BaseBuilder):
    def __init__(self, model, name: str = 'Feature'):
        super().__init__(model=model, name=name)
        self._wavelength = np.max(model.bounding_box.length)
        self._amplitude = np.min(model.bounding_box.length)
        self._centre = model.bounding_box
        self._feature = LambdaGeologicalFeature(
            function=self._function, model=self.model, name=self.name, builder=self
        )

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def wavelength(self):
        return self._wavelength

    def _function(self, xyz):
        return xyz[:, 2] + np.sin(xyz[:, 0] / self.wavelength) * self.amplitude

    def build(self, **kwargs):
        # the feature object identity is kept stable across rebuilds --
        # only its underlying function needs refreshing here.
        self._feature.function = self._function
        self._up_to_date = True

    def up_to_date(self, callback=None):
        if not self._up_to_date:
            self.update()
            if callable(callback):
                callback(1)
            return
        if callable(callback):
            callback(1)
