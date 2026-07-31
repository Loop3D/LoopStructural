"""
Finite difference interpolator with fold constraints.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from loop_common.logging import get_logger as getLogger

from ._finite_difference_interpolator import FiniteDifferenceInterpolator
from ._fold_norm_alignment import resolve_fold_norm_target
from ._fold_setup import setup_with_fold_constraints
from ._interpolatortype import InterpolatorType
from ._regularisation import DirectionalRegularisation

logger = getLogger(__name__)


_DEFAULT_FOLD_REGULARISATION = object()


from ._fold_event import FoldEvent


class FDFoldInterpolator(FiniteDifferenceInterpolator):
    """
    Finite difference interpolator that supports anisotropic fold
    regularisation.

    Analogous to :class:`DiscreteFoldInterpolator` (which works on a
    tetrahedral P1 mesh) but applied to a regular Cartesian grid.  The fold
    geometry is evaluated at every grid node and used to build direction-
    dependent second-derivative regularisation constraints via
    :meth:`minimise_directional_gradient_change`.

    The fold constraints mirror the four families in the P1 version:

    1. **Fold orientation** – gradient orthogonal to the deformed foliation.
    2. **Fold axis**        – gradient orthogonal to the fold axis.
    3. **Fold normalisation** – unit gradient magnitude in the fold-normal
       direction.
    4. **Fold regularisation** – anisotropic smoothness: stronger along the
       fold normal and fold axis, weaker across the fold.

    Parameters
    ----------
    grid
        A structured Cartesian grid (``StructuredGrid`` / rectilinear grid).
    fold : FoldEvent, optional
        Fold geometry provider with a ``get_deformed_orientation`` method.
    data : dict, optional
        Initial constraint data forwarded to the parent interpolator.
    """

    def __init__(self, grid, fold: Optional[FoldEvent] = None, data=None):
        FiniteDifferenceInterpolator.__init__(self, grid, data=data)
        self.type = InterpolatorType.FINITE_DIFFERENCE
        self.fold = fold

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def setup_interpolator(self, **kwargs):
        """
        Set up the interpolator and add fold constraints.

        Accepted keyword arguments (all optional):

        fold_weights : dict
            Per-constraint weight overrides forwarded to
            :meth:`add_fold_constraints`.  Keys are the same as the
            parameters of that method (e.g. ``fold_orientation``,
            ``fold_regularisation``, …).
        All other kwargs are forwarded to the parent
        :meth:`FiniteDifferenceInterpolator.setup_interpolator`.
        """
        return setup_with_fold_constraints(
            fold=self.fold,
            kwargs=kwargs,
            interpolator_name=self.__class__.__name__,
            base_setup=super().setup_interpolator,
            add_fold_constraints=self.add_fold_constraints,
            finalize_report=self.finalize_setup_diagnostics_report,
        )

    # ------------------------------------------------------------------
    # Fold constraints
    # ------------------------------------------------------------------

    def add_fold_constraints(
        self,
        fold_orientation: Optional[float] = 10.0,
        fold_axis_w: Optional[float] = 10.0,
        fold_regularisation=_DEFAULT_FOLD_REGULARISATION,
        fold_normalisation: Optional[float] = 1.0,
        fold_norm: Optional[float] = -1.0,
        dgz_alignment: str = "warn",
        mask_fn: Optional[Callable] = None,
    ):
        """
        Add fold geometry constraints to the finite difference system.

        Fold direction vectors are evaluated at every grid *node* by calling
        ``self.fold.get_deformed_orientation(self.support.nodes)``.  The
        returned vectors are then used to build:

        * Gradient-orientation constraints (dot product = 0) for the
          deformed foliation plane normal and the fold axis.
        * A norm constraint forcing the gradient magnitude in the fold-normal
          direction to be ``fold_norm``.
        * Anisotropic regularisation using directional second derivatives,
          applied with different weights for the three fold-geometry
          directions.

        Parameters
        ----------
        fold_orientation : float or None
            Weight for the constraint that :math:`\\nabla f` is perpendicular
            to the deformed foliation (i.e. lies in the axial plane).
            Pass ``None`` to skip.
        fold_axis_w : float or None
            Weight for the constraint that :math:`\\nabla f` is perpendicular
            to the fold axis.  Pass ``None`` to skip.
        fold_regularisation : list of three floats or None
            Weights ``[w0, w1, w2]`` for the anisotropic regularisation
            applied in the fold-normal, deformed-orientation, and fold-axis
            directions respectively.  A larger ``w0`` (fold-normal direction)
            keeps the gradient smooth *across* the fold; smaller ``w1``/``w2``
            allow it to vary more freely along the fold hinge.  Defaults to
            ``[0.1, 0.01, 0.01]``.  Pass ``None`` to skip entirely.
        fold_normalisation : float or None
            Weight for the norm constraint.  Pass ``None`` to skip.
        fold_norm : float or None
            Target gradient projection in the fold-normal direction. The
            default is ``-1.0`` to align with the common convention where
            ``dgz`` is opposite to foliation normals.
        dgz_alignment : {"none", "warn", "correct"}
            How to handle detected sign mismatch between ``dgz`` and normal
            constraints (if present). ``warn`` logs only, ``correct`` flips
            the effective ``fold_norm`` sign, and ``none`` skips checks.
        mask_fn : callable or None
            If provided, called with ``(n_nodes, 3)`` node positions; nodes
            for which the function returns ``True`` are excluded from all fold
            constraints.
        """
        if fold_regularisation is _DEFAULT_FOLD_REGULARISATION:
            fold_regularisation = [0.1, 0.01, 0.01]

        # Evaluate fold geometry at every grid node.
        node_positions = self.support.nodes  # (n_nodes, 3)
        deformed_orientation, fold_axis, dgz = self.fold.get_deformed_orientation(node_positions)
        # deformed_orientation : (n_nodes, 3) – vector lying in the axial plane
        # fold_axis            : (n_nodes, 3) – fold hinge direction
        # dgz                  : (n_nodes, 3) – fold normal (across-fold direction)

        # Build a boolean weight mask (1 = active, 0 = excluded).
        weight = np.ones(self.support.n_nodes, dtype=float)
        if mask_fn is not None:
            weight[mask_fn(node_positions)] = 0.0

        # Convenience: points array in the format expected by
        # add_gradient_orthogonal_constraints (xyz | vx vy vz | … ).
        active = weight > 0.0
        active_pos = node_positions[active]  # (n_active, 3)

        # 1. Fold orientation: ∇f · deformed_orientation = 0
        if fold_orientation is not None:
            logger.info(f"Adding fold orientation constraint  w = {fold_orientation}")
            self.add_gradient_orthogonal_constraints(
                active_pos,
                deformed_orientation[active],
                w=fold_orientation,
                name="fold orientation",
            )

        # 2. Fold axis: ∇f · fold_axis = 0
        if fold_axis_w is not None:
            logger.info(f"Adding fold axis constraint  w = {fold_axis_w}")
            self.add_gradient_orthogonal_constraints(
                active_pos,
                fold_axis[active],
                w=fold_axis_w,
                name="fold axis",
            )

        # 3. Fold normalisation: ∇f · dgz = fold_norm
        if fold_normalisation is not None:
            target_norm = resolve_fold_norm_target(
                fold=self.fold,
                normal_constraints=self.get_norm_constraints(),
                fold_norm=fold_norm,
                dgz_alignment=dgz_alignment,
                logger=logger,
                default_fold_norm=-1.0,
            )
            logger.info(f"Adding fold normalisation constraint  w = {fold_normalisation}")
            self.add_gradient_orthogonal_constraints(
                active_pos,
                dgz[active],
                w=fold_normalisation,
                b=target_norm,
                name="fold normalisation",
            )

        # 4. Anisotropic regularisation: directional second derivatives.
        if fold_regularisation is not None:
            logger.info(
                f"Adding fold regularisation  w = {fold_regularisation[0]}, "
                f"{fold_regularisation[1]}, {fold_regularisation[2]}"
            )

            def _masked_direction(vectors):
                masked = np.asarray(vectors, dtype=float).copy()
                masked[~active] = 0.0
                return masked

            self.add_directional_regularisation(
                (
                    DirectionalRegularisation(
                        weight=fold_regularisation[0],
                        direction=_masked_direction(dgz),
                        name="fold regularisation 1",
                    ),
                    DirectionalRegularisation(
                        weight=fold_regularisation[1],
                        direction=_masked_direction(deformed_orientation),
                        name="fold regularisation 2",
                    ),
                    DirectionalRegularisation(
                        weight=fold_regularisation[2],
                        direction=_masked_direction(fold_axis),
                        name="fold regularisation 3",
                    ),
                )
            )
