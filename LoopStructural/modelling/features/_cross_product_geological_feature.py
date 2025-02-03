""" """

import numpy as np
from typing import Optional

from ...modelling.features import BaseFeature

from ...utils import getLogger

logger = getLogger(__name__)


class CrossProductGeologicalFeature(BaseFeature):
    def __init__(
        self,
        name: str,
        geological_feature_a: BaseFeature,
        geological_feature_b: BaseFeature,
    ):
        """

        Create a geological feature for a vector field using the cross
        product between
        two existing features
        Parameters
        ----------
        name: feature name
        geological_feature_a: first feature
        geological_feature_b: second feature


        Parameters
        ----------
        name : str
            name of the feature
        geological_feature_a : BaseFeature
            Left hand side of cross product
        geological_feature_b : BaseFeature
            Right hand side of cross product
        """
        super().__init__(name)
        self.geological_feature_a = geological_feature_a
        self.geological_feature_b = geological_feature_b
        self.value_feature = None

    def evaluate_gradient(self, locations: np.ndarray, ignore_regions=False) -> np.ndarray:
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
        v1 = self.geological_feature_a.evaluate_gradient(locations, ignore_regions)
        # v1 /= np.linalg.norm(v1,axis=1)[:,None]
        v2 = self.geological_feature_b.evaluate_gradient(locations, ignore_regions)
        # v2 /= np.linalg.norm(v2,axis=1)[:,None]
        return np.cross(v1, v2, axisa=1, axisb=1)

    def evaluate_value(self, evaluation_points: np.ndarray, ignore_regions=False) -> np.ndarray:
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
            values[:] = self.value_feature.evaluate_value(evaluation_points, ignore_regions)
        return values

    def mean(self):
        if self.value_feature:
            return self.value_feature.mean()
        return 0.0

    def min(self):
        if self.value_feature:
            return self.value_feature.min()
        return 0.0

    def max(self):
        if self.value_feature:
            return self.value_feature.max()
        return 0.0

    def get_data(self, value_map: Optional[dict] = None):
        return

    def copy(self, name: Optional[str] = None):
        if name is None:
            name = f'{self.name}_copy'
        return CrossProductGeologicalFeature(
            name, self.geological_feature_a, self.geological_feature_b
        )
