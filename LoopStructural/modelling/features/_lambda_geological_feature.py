"""
Geological features
"""
from LoopStructural.utils.maths import regular_tetraherdron_for_points, gradient_from_tetrahedron
from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
import numpy as np
from typing import Callable, Optional
from ...utils import LoopValueError

logger = getLogger(__name__)


class LambdaGeologicalFeature(BaseFeature):
    def __init__(
        self,
        function: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        name: str = "unnamed_lambda",
        gradient_function: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        model=None,
        regions: Optional[list] = None,
        faults: Optional[list] = None,
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
        BaseFeature.__init__(self, name, model, faults if faults is not None else [], regions if regions is not None else [], builder)
        self.type = FeatureType.LAMBDA
        self.function = function
        self.gradient_function = gradient_function
        self.regions = regions if regions is not None else []

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
        v[:] = np.nan

        mask = self._calculate_mask(pos, ignore_regions=ignore_regions)
        pos = self._apply_faults(pos)
        if self.function is not None:
            v[mask] = self.function(pos[mask,:])
        return v

    def evaluate_gradient(self, pos: np.ndarray, ignore_regions=False,element_scale_parameter=None) -> np.ndarray:
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
        if pos.shape[1] != 3:
            raise LoopValueError("Need Nx3 array of xyz points to evaluate gradient")
        logger.info(f'Calculating gradient for {self.name}')
        if element_scale_parameter is None:
            if self.model is not None:
                element_scale_parameter = np.min(self.model.bounding_box.step_vector) / 10
            else:
                element_scale_parameter = 1
        else:
            try:
                element_scale_parameter = float(element_scale_parameter)
            except ValueError:
                logger.error("element_scale_parameter must be a float")
                element_scale_parameter = 1
        v = np.zeros((pos.shape[0], 3))
        v = np.zeros(pos.shape)
        v[:] = np.nan
        mask = self._calculate_mask(pos, ignore_regions=ignore_regions)
        # evaluate the faults on the nodes of the faulted feature support
        # then evaluate the gradient at these points
        if len(self.faults) > 0:
            # generate a regular tetrahedron for each point
            # we will then move these points by the fault and then recalculate the gradient.
            # this should work...
            resolved = False
            tetrahedron = regular_tetraherdron_for_points(pos, element_scale_parameter)

            while resolved:
                for f in self.faults:
                    v = (
                        f[0]
                        .evaluate_value(tetrahedron.reshape(-1, 3), fillnan='nearest')
                        .reshape(tetrahedron.shape[0], 4)
                    )
                    flag = np.logical_or(np.all(v > 0, axis=1), np.all(v < 0, axis=1))
                    if np.any(~flag):
                        logger.warning(
                            f"Points are too close to fault {f[0].name}. Refining the tetrahedron"
                        )
                        element_scale_parameter *= 0.5
                        tetrahedron = regular_tetraherdron_for_points(pos, element_scale_parameter)

                resolved = True

            tetrahedron_faulted = self._apply_faults(np.array(tetrahedron.reshape(-1, 3))).reshape(
                tetrahedron.shape
            )

            values = self.function(tetrahedron_faulted.reshape(-1, 3)).reshape(
                (-1, 4)
            )
            v[mask, :] = gradient_from_tetrahedron(tetrahedron[mask, :, :], values[mask])

            return v
        if self.gradient_function is None:
            v[:, :] = np.nan
        else:
            v[:, :] = self.gradient_function(pos)
        return v

    def get_data(self, value_map: Optional[dict] = None):
        return

    def copy(self, name: Optional[str] = None):
        return LambdaGeologicalFeature(
            self.function,
            name if name is not None else f'{self.name}_copy',
            self.gradient_function,
            self.model,
            self.regions,
            self.faults,
            self.builder,
        )
    def is_valid(self):
        if self.function is None and self.gradient_function is None:
            return False
        return True
