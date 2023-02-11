"""
Geological features
"""
from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
import numpy as np

logger = getLogger(__name__)


class LambdaGeologicalFeature(BaseFeature):
    def __init__(
        self,
        function=None,
        name="unnamed_lambda",
        gradient_function=None,
        model=None,
        regions=[],
        faults=[],
        builder=None,
    ):
        """_summary_

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

    def evaluate_value(self, xyz: np.ndarray) -> np.ndarray:
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
        v = np.zeros((xyz.shape[0]))
        if self.function is None:
            v[:] = np.nan
        else:
            v[:] = self.function(xyz)
        return v

    def evaluate_gradient(self, xyz: np.ndarray) -> np.ndarray:
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
        v = np.zeros((xyz.shape[0], 3))
        if self.gradient_function is None:
            v[:, :] = np.nan
        else:
            v[:, :] = self.gradient_function(xyz)
        return v
