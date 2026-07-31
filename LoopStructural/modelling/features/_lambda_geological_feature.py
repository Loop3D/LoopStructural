"""
Geological features
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from LoopStructural.utils.maths import gradient_from_tetrahedron, regular_tetraherdron_for_points

from ...modelling.features import BaseFeature, FeatureType
from ...utils import LoopValueError, getLogger

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
        function : Callable[[np.ndarray], np.ndarray], optional
            function that takes an Nx3 array of xyz points and returns the value of the
            feature at each point, by default None
        name : str, optional
            name of the feature, by default "unnamed_lambda"
        gradient_function : Callable[[np.ndarray], np.ndarray], optional
            function that takes an Nx3 array of xyz points and returns the gradient of the
            feature at each point, by default None
        model : GeologicalModel, optional
            the geological model this feature is associated with, by default None
        regions : list, optional
            list of regions to restrict where this feature is evaluated, by default []
        faults : list, optional
            list of faults that affect this feature, by default []
        builder : optional
            the builder used to create this feature, by default None
        """
        BaseFeature.__init__(self, name, model, faults if faults is not None else [], regions if regions is not None else [], builder)
        self.type = FeatureType.LAMBDA
        self.function = function
        self.gradient_function = gradient_function
        self.regions = regions if regions is not None else []

    def evaluate_value(self, pos: np.ndarray, ignore_regions=False) -> np.ndarray:
        """Evaluate the value of the underlying function at locations, applying
        any faults and regions associated with this feature

        Parameters
        ----------
        pos : np.ndarray
            Nx3 array of xyz locations to evaluate the feature at
        ignore_regions : bool, optional
            whether to ignore the regions associated with this feature, by default False

        Returns
        -------
        np.ndarray
            value of the feature at each location, nan where outside of the regions
        """
        v = np.zeros(pos.shape[0])
        v[:] = np.nan

        # Precompute each fault's scalar value (gx = fault.__getitem__(0).evaluate_value)
        # once and reuse for both the region mask check and fault application.
        # FaultSegment.evaluate_value(pos) == self.__getitem__(0).evaluate_value(pos) == gx.
        # apply_to_points also needs gx to determine the displacement mask — passing it
        # avoids the duplicate evaluation on the np.copy of pos.
        precomputed_gx = {id(f): f.evaluate_value(pos) for f in self.faults}

        # Region mask: pass precomputed gx so PositiveRegion/NegativeRegion skip re-evaluation
        mask = np.ones(pos.shape[0], dtype=bool)
        if not ignore_regions:
            for r in self.regions:
                pre = precomputed_gx.get(id(getattr(r, 'feature', None)))
                mask = np.logical_and(mask, r(pos, precomputed_val=pre))

        # Apply faults: pass precomputed gx so apply_to_points skips one evaluate_value call
        if self.faults_enabled:
            for f in self.faults:
                pos = f.apply_to_points(pos, precomputed_gx=precomputed_gx.get(id(f)))

        if self.function is not None:
            v[mask] = self.function(pos[mask,:])
        return v

    def evaluate_gradient(self, pos: np.ndarray, ignore_regions=False,element_scale_parameter=None) -> np.ndarray:
        """Evaluate the gradient of the underlying function at locations, applying
        any faults associated with this feature

        Parameters
        ----------
        pos : np.ndarray
            Nx3 array of xyz locations to evaluate the gradient at
        ignore_regions : bool, optional
            whether to ignore the regions associated with this feature, by default False
        element_scale_parameter : float, optional
            size of the finite tetrahedron used to numerically estimate the gradient when
            faults are present, by default a tenth of the model's minimum step vector

        Returns
        -------
        np.ndarray
            Nx3 array of the gradient of the feature at each location, nan where undefined
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
