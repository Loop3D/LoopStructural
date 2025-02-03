""" """

import numpy as np
from typing import Optional

from ...modelling.features import BaseFeature

from ...utils import getLogger

logger = getLogger(__name__)


class ProjectedVectorFeature(BaseFeature):
    def __init__(
        self,
        name: str,
        vector: np.ndarray,
        plane_feature: BaseFeature,
    ):
        """

        Create a geological feature by projecting a vector onto a feature representing a plane
        E.g. project a thickness vector onto an axial surface

        Parameters
        ----------
        name: feature name
        vector: the vector to project
        plane_feature: the plane


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
        self.plane_feature = plane_feature
        self.vector = vector

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

        # project s0 onto axis plane B X A X B
        plane_normal = self.plane_feature.evaluate_gradient(locations, ignore_regions)
        vector = np.tile(self.vector, (locations.shape[0], 1))

        projected_vector = np.cross(
            plane_normal, np.cross(vector, plane_normal, axisa=1, axisb=1), axisa=1, axisb=1
        )
        return projected_vector

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
        return ProjectedVectorFeature(
            name=name, vector=self.vector, plane_feature=self.plane_feature
        )
