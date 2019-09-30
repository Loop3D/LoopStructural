import numpy as np

from LoopStructural.modelling.features.geological_feature import GeologicalFeature


class CompositeGeologicalFeature(GeologicalFeature):
    def __init__(self, geological_feature_a, region_a, geological_feature_b, region_b):
        self.geological_feature_a = geological_feature_a
        self.region_a = region_a

        self.geological_feature_b = geological_feature_b
        self.region_b = region_b

    def evaluate_gradient(self, locations):

        pass

    def evaluate_value(self, evaluation_points):

        pass

    def mean(self):
        return np.nanmean(self.support.get_node_values())

    def min(self):
        return np.nanmin(self.support.get_node_values())

    def max(self):
        return np.nanmax(self.support.get_node_values())