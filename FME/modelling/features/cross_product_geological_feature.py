import numpy as np

from FME.modelling.features.geological_feature import GeologicalFeature


class CrossProductGeologicalFeature(GeologicalFeature):
    def __init__(self, name, geological_feature_a, geological_feature_b):
        super().__init__(name+"_faulted", geological_feature_a.support)
        self.geological_feature_a = geological_feature_a
        self.geological_feature_b = geological_feature_b
    def evaluate_gradient(self, locations):
        return np.cross(self.geological_feature_a.evaluate_gradient(locations),
                        self.geological_feature_b.evaluate_gradient(locations),
                        axisa=1,
                        axisb=1)
    def evaluate_value(self, evaluation_points):
        return 0.

    def mean(self):
        return 0.

    def min(self):
        return 0.

    def max(self):
        return 0.