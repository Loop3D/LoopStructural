"""
"""
import logging

import numpy as np

from LoopStructural.modelling.features.geological_feature import \
    GeologicalFeature

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class CrossProductGeologicalFeature(GeologicalFeature):
    """[summary]

    [extended_summary]

    Parameters
    ----------
    GeologicalFeature : [type]
        [description]
    """
    def __init__(self, name, geological_feature_a, geological_feature_b):
        """
        Create a geological feature for a vector field using the cross
        product between
        two existing features
        Parameters
        ----------
        name: feature name
        geological_feature_a: first feature
        geological_feature_b: second feature
        """
        super().__init__(name, None)
        self.geological_feature_a = geological_feature_a
        self.geological_feature_b = geological_feature_b
        self.value_feature = None
    def evaluate_gradient(self, locations):
        """
        Calculate the gradient of the geological feature by using numpy to
        calculate the cross
        product between the two existing feature gradients.
        This means both features have to be evaluated for the locations
        Parameters
        ----------
        locations

        Returns
        -------

        """
        v1 = self.geological_feature_a.evaluate_gradient(locations)
        # v1 /= np.linalg.norm(v1,axis=1)[:,None]
        v2 = self.geological_feature_b.evaluate_gradient(locations)
        # v2 /= np.linalg.norm(v2,axis=1)[:,None]
        return np.cross(v1,v2,
                        axisa=1,
                        axisb=1)

    def evaluate_value(self, evaluation_points):
        """
        Return 0 because there is no value for this feature
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        values = np.zeros(evaluation_points.shape[0])
        if self.value_feature is not None:
            values[:] = self.value_feature.evaluate_value(evaluation_points)
        return values

    def mean(self):
        if self.value_feature:
            return self.value_feature.mean()
        return 0.

    def min(self):
        if self.value_feature:
            return self.value_feature.min()
        return 0.

    def max(self):
        if self.value_feature:
            return self.value_feature.max()
        return 0.
