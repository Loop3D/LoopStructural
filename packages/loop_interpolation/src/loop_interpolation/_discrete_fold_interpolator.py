"""
Piecewise linear interpolator using folds
"""

from typing import Optional, Callable

import numpy as np

from ._p1interpolator import P1Interpolator as PiecewiseLinearInterpolator
from ._interpolatortype import InterpolatorType
from ._regularisation import DirectionalRegularisation
from ._fold_setup import setup_with_fold_constraints
from ._fold_norm_alignment import resolve_fold_norm_target

from loop_common.logging import get_logger as getLogger
from loop_common.math import rng
from ._fold_event import FoldEvent  # noqa: F401  (re-exported for convenience)

logger = getLogger(__name__)


class DiscreteFoldInterpolator(PiecewiseLinearInterpolator):
    """ """

    def __init__(self, support, fold: Optional[FoldEvent] = None):
        """
        A piecewise linear interpolator that can also use fold constraints defined in Laurent et al., 2016

        Parameters
        ----------
        support
            discrete support with nodes and elements etc
        fold FoldEvent
            a fold event with a valid geometry
        """

        PiecewiseLinearInterpolator.__init__(self, support)
        self.type = InterpolatorType.DISCRETE_FOLD
        self.fold = fold

    def update_fold(self, fold):
        """

        Parameters
        ----------
        fold : FoldEvent
            a fold that contrains the geometry we are trying to add

        Returns
        -------

        """
        logger.error("updating fold, this should be done by accessing the fold attribute")
        self.fold = fold

    def setup_interpolator(self, **kwargs):
        return setup_with_fold_constraints(
            fold=self.fold,
            kwargs=kwargs,
            interpolator_name=self.__class__.__name__,
            base_setup=super().setup_interpolator,
            add_fold_constraints=self.add_fold_constraints,
            finalize_report=self.finalize_setup_diagnostics_report,
        )

    def add_fold_constraints(
        self,
        fold_orientation=10.0,
        fold_axis_w=10.0,
        fold_regularisation=(0.1, 0.01, 0.01),
        fold_normalisation=1.0,
        fold_norm=-1.0,
        dgz_alignment="warn",
        step=2,
        mask_fn: Optional[Callable] = None,
    ):
        """

        Parameters
        ----------
        fold_orientation : double
            weight for the fold direction/orientation in the least squares system
        fold_axis_w : double
            weight for the fold axis in the least squares system
        fold_regularisation : list
            weight for the fold regularisation in the least squares system
        fold_normalisation : double
            weight for the fold norm constraint in the least squares system
        fold_norm
            length of the interpolation norm in the least squares system
        dgz_alignment : {"none", "warn", "correct"}
            How to handle detected sign mismatch between ``dgz`` and normal
            constraints (if present). ``warn`` logs only, ``correct`` flips
            the effective ``fold_norm`` sign, and ``none`` skips checks.
        step: int
            array step for adding constraints


        Returns
        -------

        Notes
        -----
        For more information about the fold weights see EPSL paper by Gautier Laurent 2016

        """
        # get the gradient of all of the elements of the mesh
        eg = self.support.get_element_gradients(np.arange(self.support.n_elements))
        # get array of all nodes for all elements N,4,3
        nodes = self.support.nodes[self.support.get_elements()[np.arange(self.support.n_elements)]]
        # calculate the fold geometry for the elements barycentre
        deformed_orientation, fold_axis, dgz = self.fold.get_deformed_orientation(
            self.support.barycentre
        )
        element_idx = np.arange(self.support.n_elements)
        rng.shuffle(element_idx)
        # calculate element volume for weighting
        vecs = nodes[:, 1:, :] - nodes[:, 0, None, :]
        vol = np.abs(np.linalg.det(vecs)) / 6
        weight = np.ones(self.support.n_elements, dtype=float)
        if mask_fn is not None:
            weight[mask_fn(self.support.barycentre)] = 0
        if fold_orientation is not None:
            """
            dot product between vector in deformed ori plane = 0
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold orientation constraint to w = {fold_orientation}")
            selected_idx = element_idx[::step]
            A = np.einsum(
                "ij,ijk->ik",
                deformed_orientation[selected_idx, :],
                eg[selected_idx, :, :],
            )
            A *= vol[selected_idx, None]
            B = np.zeros(A.shape[0])
            idc = self.support.get_elements()[selected_idx, :]
            self.add_constraints_to_least_squares(
                A,
                B,
                idc,
                w=weight[selected_idx] * fold_orientation,
                name="fold orientation",
            )

        if fold_axis_w is not None:
            """
            dot product between axis and gradient should be 0
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold axis constraint to  w = {fold_axis_w}")
            selected_idx = element_idx[::step]
            A = np.einsum(
                "ij,ijk->ik",
                fold_axis[selected_idx, :],
                eg[selected_idx, :, :],
            )
            A *= vol[selected_idx, None]
            B = np.zeros(A.shape[0]).tolist()
            idc = self.support.get_elements()[selected_idx, :]

            self.add_constraints_to_least_squares(
                A,
                B,
                idc,
                w=weight[selected_idx] * fold_axis_w,
                name="fold axis",
            )

        if fold_normalisation is not None:
            """
            specify scalar norm in X direction
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold normalisation constraint to  w = {fold_normalisation}")
            selected_idx = element_idx[::step]
            A = np.einsum("ij,ijk->ik", dgz[selected_idx, :], eg[selected_idx, :, :])
            A *= vol[selected_idx, None]

            target_norm = resolve_fold_norm_target(
                fold=self.fold,
                normal_constraints=self.get_norm_constraints(),
                fold_norm=fold_norm,
                dgz_alignment=dgz_alignment,
                logger=logger,
                default_fold_norm=-1.0,
            )
            B = np.ones(A.shape[0]) * target_norm
            B *= fold_normalisation
            B *= vol[selected_idx]
            idc = self.support.get_elements()[selected_idx, :]

            self.add_constraints_to_least_squares(
                A,
                B,
                idc,
                w=weight[selected_idx] * fold_normalisation,
                name="fold normalisation",
            )

        if fold_regularisation is not None:
            """
            fold constant gradient
            """
            logger.info(
                f"Adding fold regularisation constraint to  w = {fold_regularisation[0]} {fold_regularisation[1]} {fold_regularisation[2]}"
            )

            def _masked_fold_direction(component_index: int):
                def _provider(points: np.ndarray) -> np.ndarray:
                    deformed, axis, normal = self.fold.get_deformed_orientation(points)
                    vectors = (normal, deformed, axis)[component_index]
                    if mask_fn is not None:
                        masked = np.asarray(vectors, dtype=float).copy()
                        masked[mask_fn(points)] = 0.0
                        return masked
                    return vectors

                return _provider

            self.add_directional_regularisation(
                (
                    DirectionalRegularisation(
                        weight=fold_regularisation[0],
                        direction=_masked_fold_direction(0),
                        name="fold regularisation 1",
                    ),
                    DirectionalRegularisation(
                        weight=fold_regularisation[1],
                        direction=_masked_fold_direction(1),
                        name="fold regularisation 2",
                    ),
                    DirectionalRegularisation(
                        weight=fold_regularisation[2],
                        direction=_masked_fold_direction(2),
                        name="fold regularisation 3",
                    ),
                )
            )
