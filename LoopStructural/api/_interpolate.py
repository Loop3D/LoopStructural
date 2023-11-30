from ast import LShift
import numpy as np
import pandas as pd

from typing import Optional
from LoopStructural.interpolators import (
    GeologicalInterpolator,
    InterpolatorFactory,
)
from LoopStructural.utils import BoundingBox
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class LoopInterpolator:
    def __init__(
        self,
        bounding_box: BoundingBox,
        dimensions: int = 3,
        type: str = "FDI",
        nelements: int = 1000,
    ):
        """Scikitlearn like interface for LoopStructural interpolators
        useful for quickly building an interpolator to apply to a dataset
        build a generic interpolation object speficying the bounding box
        and then fit to constraints and evaluate the interpolator

        Parameters
        ----------
        bounding_box : BoundingBox
            Area where the interpolation will work
        dimensions : int, optional
            number of dimensions e.g. 3d or 2d, by default 3
        type : str, optional
            type of interpolation algorithm FDI- finite difference, PLI - linear finite elements,
              by default "FDI"
        nelements : int, optional
            degrees of freedom for interpolator, by default 1000
        """
        self.dimensions = dimensions
        self.type = "FDI"
        self.interpolator: GeologicalInterpolator = (
            InterpolatorFactory.create_interpolator(
                type,
                bounding_box,
                nelements,
            )
        )

    def fit(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_constraints: Optional[np.ndarray] = None,
    ):
        """_summary_

        Parameters
        ----------
        values : Optional[np.ndarray], optional
            Value constraints for implicit function, by default None
        tangent_vectors : Optional[np.ndarray], optional
            tangent constraints for implicit function, by default None
        normal_vectors : Optional[np.ndarray], optional
            gradient norm constraints for implicit function, by default None
        inequality_constraints : Optional[np.ndarray], optional
            _description_, by default None
        """
        if values is not None:
            self.interpolator.set_value_constraints(values)
        if tangent_vectors is not None:
            self.interpolator.set_tangent_constraints(tangent_vectors)
        if normal_vectors is not None:
            self.interpolator.set_normal_constraints(normal_vectors)
        if inequality_constraints:
            pass

        self.interpolator.setup()

    def evaluate_scalar_value(self, locations: np.ndarray) -> np.ndarray:
        """Evaluate the value of the interpolator at locations

        Parameters
        ----------
        locations : np.ndarray
            Nx3 array of locations to evaluate the interpolator at

        Returns
        -------
        np.ndarray
            value of implicit function at locations
        """
        self.interpolator.update()
        return self.interpolator.evaluate_value(locations)

    def evaluate_gradient(self, locations: np.ndarray) -> np.ndarray:
        """Evaluate the gradient of the interpolator at locations

        Parameters
        ----------
        locations : np.ndarray
            Nx3 locations

        Returns
        -------
        np.ndarray
            Nx3 gradient of implicit function
        """
        self.interpolator.update()
        return self.interpolator.evaluate_gradient(locations)

    def fit_and_evaluate_value(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_constraints: Optional[np.ndarray] = None,
    ):
        # get locations
        self.fit(
            values=values,
            tangent_vectors=tangent_vectors,
            normal_vectors=normal_vectors,
            inequality_constraints=inequality_constraints,
        )
        locations = self.interpolator.get_data_locations()
        return self.evalute_scalar_value(locations)

    def fit_and_evaluate_gradient(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_constraints: Optional[np.ndarray] = None,
    ):
        self.fit(
            values=values,
            tangent_vectors=tangent_vectors,
            normal_vectors=normal_vectors,
            inequality_constraints=inequality_constraints,
        )
        locations = self.interpolator.get_data_locations()
        return self.evaluate_gradient(locations)
