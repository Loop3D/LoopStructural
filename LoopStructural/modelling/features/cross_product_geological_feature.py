"""
"""
import logging

import numpy as np

from LoopStructural.modelling.features.geological_feature import \
    GeologicalFeature

logger = logging.getLogger(__name__)


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
        return np.cross(self.geological_feature_a.evaluate_gradient(locations),
                        self.geological_feature_b.evaluate_gradient(locations),
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
        return np.zeros(evaluation_points.shape[0])

    def mean(self):
        return 0.

    def min(self):
        return 0.

    def max(self):
        return 0.
