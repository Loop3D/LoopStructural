import numpy as np

from typing import Optional
from LoopStructural.interpolators import (
    GeologicalInterpolator,
    InterpolatorFactory,
    InterpolatorType,
)
from LoopStructural.datatypes import BoundingBox
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class LoopInterpolator:
    def __init__(
        self,
        bounding_box: BoundingBox,
        dimensions: int = 3,
        type=InterpolatorType.FINITE_DIFFERENCE,
        nelements: int = 1000,
        interpolator_setup_kwargs={},
        buffer: float = 0.2,
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
        logger.warning("LoopInterpolator is experimental and the API is subject to change")
        self.dimensions = dimensions
        self.type = "FDI"
        self.bounding_box = bounding_box
        self.interpolator: GeologicalInterpolator = InterpolatorFactory.create_interpolator(
            type, bounding_box, nelements, buffer=buffer
        )
        self.interpolator_setup_kwargs = interpolator_setup_kwargs

    def fit(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_value_constraints: Optional[np.ndarray] = None,
        inequality_pairs_constraints: Optional[np.ndarray] = None,
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
        if inequality_value_constraints is not None:
            self.interpolator.set_value_inequality_constraints(inequality_value_constraints)
        if inequality_pairs_constraints is not None:
            self.interpolator.set_inequality_pairs_constraints(inequality_pairs_constraints)
        self.interpolator.setup(**self.interpolator_setup_kwargs)

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
        inequality_value_constraints: Optional[np.ndarray] = None,
        inequality_pairs_constraints: Optional[np.ndarray] = None,
    ):
        # get locations
        self.fit(
            values=values,
            tangent_vectors=tangent_vectors,
            normal_vectors=normal_vectors,
            inequality_value_constraints=inequality_value_constraints,
            inequality_pairs_constraints=inequality_pairs_constraints,
        )
        locations = self.interpolator.get_data_locations()
        return self.evaluate_scalar_value(locations)

    def fit_and_evaluate_gradient(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_value_constraints: Optional[np.ndarray] = None,
        inequality_pairs_constraints: Optional[np.ndarray] = None,
    ):
        self.fit(
            values=values,
            tangent_vectors=tangent_vectors,
            normal_vectors=normal_vectors,
            inequality_value_constraints=inequality_value_constraints,
            inequality_pairs_constraints=inequality_pairs_constraints,
        )
        locations = self.interpolator.get_data_locations()
        return self.evaluate_gradient(locations)

    def fit_and_evaluate_value_and_gradient(
        self,
        values: Optional[np.ndarray] = None,
        tangent_vectors: Optional[np.ndarray] = None,
        normal_vectors: Optional[np.ndarray] = None,
        inequality_value_constraints: Optional[np.ndarray] = None,
        inequality_pairs_constraints: Optional[np.ndarray] = None,
    ):
        self.fit(
            values=values,
            tangent_vectors=tangent_vectors,
            normal_vectors=normal_vectors,
            inequality_value_constraints=inequality_value_constraints,
            inequality_pairs_constraints=inequality_pairs_constraints,
        )
        locations = self.interpolator.get_data_locations()
        return self.evaluate_scalar_value(locations), self.evaluate_gradient(locations)

    def plot(self, ax=None, **kwargs):
        """Plots a 2d map scalar field or 3d pyvista plot

        Parameters
        ----------
        ax : matplotlib axes, optional
            The axes you want to add the plot to, otherwise it will make one, by default None

        Returns
        -------
        _type_
            _description_
        """
        if self.dimensions == 3:
            vtkgrid = self.interpolator.support.vtk()
            vtkgrid['val'] = self.interpolator.c
            vtkgrid.plot(**kwargs)
            return vtkgrid
        elif self.dimensions == 2:
            if ax is None:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots()
            val = self.interpolator.c
            val = np.rot90(val.reshape(self.interpolator.support.nsteps, order='F'), 3)
            ax.imshow(
                val,
                origin='lower',
                extent=[
                    self.bounding_box.origin[0],
                    self.bounding_box.maximum[0],
                    self.bounding_box.origin[1],
                    self.bounding_box.maximum[1],
                ],
                **kwargs,
            )
            return val, ax

    # def isovalue(self, value: float, **kwargs):
    #     self.interpolator.isovalue(value, **kwargs)
