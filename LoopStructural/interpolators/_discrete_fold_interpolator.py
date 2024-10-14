"""
Piecewise linear interpolator using folds
"""

from typing import Optional, Callable

import numpy as np

from ..interpolators import PiecewiseLinearInterpolator, InterpolatorType
from ..modelling.features.fold import FoldEvent
from ..utils import getLogger, rng

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
        if self.fold is None:
            raise Exception("No fold event specified")
        fold_weights = kwargs.get("fold_weights", {})
        super().setup_interpolator(**kwargs)
        self.add_fold_constraints(**fold_weights)

        return

    def add_fold_constraints(
        self,
        fold_orientation=10.0,
        fold_axis_w=10.0,
        fold_regularisation=[0.1, 0.01, 0.01],
        fold_normalisation=1.0,
        fold_norm=1.0,
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
        weight = weight[::step]
        if fold_orientation is not None:
            """
            dot product between vector in deformed ori plane = 0
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold orientation constraint to w = {fold_orientation}")
            A = np.einsum(
                "ij,ijk->ik",
                deformed_orientation[element_idx[::step], :],
                eg[element_idx[::step], :, :],
            )
            A *= vol[element_idx[::step], None]
            B = np.zeros(A.shape[0])
            idc = self.support.get_elements()[element_idx[::step], :]
            self.add_constraints_to_least_squares(
                A, B, idc, w=weight * fold_orientation, name="fold orientation"
            )

        if fold_axis_w is not None:
            """
            dot product between axis and gradient should be 0
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold axis constraint to  w = {fold_axis_w}")
            A = np.einsum(
                "ij,ijk->ik",
                fold_axis[element_idx[::step], :],
                eg[element_idx[::step], :, :],
            )
            A *= vol[element_idx[::step], None]
            B = np.zeros(A.shape[0]).tolist()
            idc = self.support.get_elements()[element_idx[::step], :]

            self.add_constraints_to_least_squares(
                A, B, idc, w=weight * fold_axis_w, name="fold axis"
            )

        if fold_normalisation is not None:
            """
            specify scalar norm in X direction
            """
            rng.shuffle(element_idx)

            logger.info(f"Adding fold normalisation constraint to  w = {fold_normalisation}")
            A = np.einsum("ij,ijk->ik", dgz[element_idx[::step], :], eg[element_idx[::step], :, :])
            A *= vol[element_idx[::step], None]

            B = np.ones(A.shape[0])

            if fold_norm is not None:
                B[:] = fold_norm
            B *= fold_normalisation
            B *= vol[element_idx[::step]]
            idc = self.support.get_elements()[element_idx[::step], :]

            self.add_constraints_to_least_squares(
                A, B, idc, w=weight * fold_normalisation, name="fold normalisation"
            )

        if fold_regularisation is not None:
            """
            fold constant gradient
            """
            logger.info(
                f"Adding fold regularisation constraint to  w = {fold_regularisation[0]} {fold_regularisation[1]} {fold_regularisation[2]}"
            )
            self.minimise_edge_jumps(
                w=fold_regularisation[0], vector=dgz, name="fold regularisation 1"
            )
            self.minimise_edge_jumps(
                w=fold_regularisation[1],
                vector=deformed_orientation,
                name="fold regularisation 2",
            )
            self.minimise_edge_jumps(
                w=fold_regularisation[2], vector=fold_axis, name="fold regularisation 3"
            )
