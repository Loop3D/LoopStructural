"""
Geological features
"""

from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
import numpy as np
from typing import Callable, Optional

logger = getLogger(__name__)


class LambdaGeologicalFeature(BaseFeature):

    def __init__(
        self,
        function: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        name: str = "unnamed_lambda",
        gradient_function: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        model=None,
        regions: list = [],
        faults: list = [],
        builder=None,
    ):
        """A lambda geological feature is a wrapper for a geological
        feature that has a function at the base. This can be then used
        in place of a geological feature.

        Parameters
        ----------
        function : _type_, optional
            _description_, by default None
        name : str, optional
            _description_, by default "unnamed_lambda"
        gradient_function : _type_, optional
            _description_, by default None
        model : _type_, optional
            _description_, by default None
        regions : list, optional
            _description_, by default []
        faults : list, optional
            _description_, by default []
        builder : _type_, optional
            _description_, by default None
        """
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.type = FeatureType.LAMBDA
        self.function = function
        self.gradient_function = gradient_function

    def evaluate_value(self, pos: np.ndarray, ignore_regions=False) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        xyz : np.ndarray
            _description_

        Returns
        -------
        np.ndarray
            _description_
        """
        v = np.zeros((pos.shape[0]))
        if self.function is None:
            v[:] = np.nan
        else:
            v[:] = self.function(pos)
        return v

    def evaluate_gradient(self, pos: np.ndarray, ignore_regions=False) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        xyz : np.ndarray
            _description_

        Returns
        -------
        np.ndarray
            _description_
        """
        v = np.zeros((pos.shape[0], 3))
        if self.gradient_function is None:
            v[:, :] = np.nan
        else:
            v[:, :] = self.gradient_function(pos)
        return v

    def get_data(self, value_map: Optional[dict] = None):
        return
