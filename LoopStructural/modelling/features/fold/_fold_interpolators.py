"""Fold-aware interpolator built on loop_interpolation's generic anisotropy constraints.

The locally-varying-anisotropy mathematics (direction-dependent gradient
orthogonality/regularisation constraints) live in
:func:`loop_interpolation.add_element_anisotropy_constraints`, a plain function
built entirely out of :class:`~loop_interpolation.P1Interpolator`'s existing
generic constraint methods — it knows nothing about folds. This module
supplies the fold-specific glue: an adapter from
:class:`~LoopStructural.modelling.features.fold.FoldEvent`'s geological
vocabulary (``get_deformed_orientation``) to the generic ``get_directions``
direction-provider protocol, and a thin ``P1Interpolator`` subclass that
keeps the familiar fold vocabulary (``fold_orientation``, ``fold_axis_w``,
…) for the rest of LoopStructural.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Callable

import numpy.typing as npt
from loop_common.supports._base_support import BaseSupport
from loop_interpolation import (
    ConstraintDiagnosticsReport,
    InterpolatorType,
    P1Interpolator,
    add_element_anisotropy_constraints,
    interpolator_map,
)

from ....utils import getLogger

if TYPE_CHECKING:
    from ._fold import FoldEvent

logger = getLogger(__name__)


class _FoldDirectionAdapter:
    """Adapts a FoldEvent to the generic ``get_directions`` protocol."""

    def __init__(self, fold: FoldEvent) -> None:
        """Wrap a fold event so it satisfies the generic direction-provider protocol.

        Parameters
        ----------
        fold
            Fold event exposing ``get_deformed_orientation(points)``.
        """
        self.fold = fold

    def get_directions(
        self, points: npt.NDArray
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
        """Return the ``(primary, secondary, normal)`` direction fields at *points*."""
        return self.fold.get_deformed_orientation(points)


class DiscreteFoldInterpolator(P1Interpolator):
    """Piecewise linear interpolator that supports fold constraints."""

    def __init__(self, support: BaseSupport, fold: FoldEvent | None = None) -> None:
        """Create a piecewise linear interpolator that can use fold constraints from Laurent et al., 2016.

        Parameters
        ----------
        support
            discrete support with nodes and elements etc
        fold : FoldEvent, optional
            a fold event with a valid geometry
        """
        P1Interpolator.__init__(self, support)
        self.type = InterpolatorType.DISCRETE_FOLD
        self._fold = fold

    @property
    def fold(self) -> FoldEvent | None:
        """The fold event providing this interpolator's constraint geometry."""
        return self._fold

    @fold.setter
    def fold(self, value: FoldEvent | None) -> None:
        self._fold = value

    def update_fold(self, fold: FoldEvent | None) -> None:
        """Replace the fold event associated with this interpolator.

        Parameters
        ----------
        fold : FoldEvent
            a fold that contains the geometry we are trying to add
        """
        logger.error("updating fold, this should be done by accessing the fold attribute")
        self.fold = fold

    def setup_interpolator(self, **kwargs: object) -> ConstraintDiagnosticsReport:
        """Set up the interpolator, then add fold constraints from ``fold_weights``."""
        if self._fold is None:
            raise RuntimeError(
                "DiscreteFoldInterpolator: no fold event set. Assign self.fold before "
                "calling setup_interpolator."
            )
        setup_kwargs = dict(kwargs)
        fold_weights = setup_kwargs.pop("fold_weights", {})
        P1Interpolator.setup_interpolator(self, **setup_kwargs)
        self.add_fold_constraints(**fold_weights)
        return self.finalize_setup_diagnostics_report()

    def add_fold_constraints(
        self,
        fold_orientation: float | None = 10.0,
        fold_axis_w: float | None = 10.0,
        fold_regularisation: Sequence[float] | None = (0.1, 0.01, 0.01),
        fold_normalisation: float | None = 1.0,
        fold_norm: float | None = -1.0,
        dgz_alignment: str = "warn",
        step: int = 2,
        mask_fn: Callable | None = None,
    ) -> None:
        """Add fold geometry constraints to the least squares system.

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
        step : int
            array step for adding constraints

        Notes
        -----
        For more information about the fold weights see EPSL paper by Gautier Laurent 2016.
        """
        add_element_anisotropy_constraints(
            self,
            _FoldDirectionAdapter(self._fold),
            primary_w=fold_orientation,
            secondary_w=fold_axis_w,
            regularisation_weights=fold_regularisation,
            normal_w=fold_normalisation,
            norm_target=fold_norm,
            norm_alignment=dgz_alignment,
            step=step,
            mask_fn=mask_fn,
        )


# Register the fold-aware wrapper as the concrete class InterpolatorFactory
# builds for InterpolatorType.DISCRETE_FOLD ("DFI"), in place of the plain
# generic P1Interpolator. interpolator_map is the same dict object
# loop_interpolation.InterpolatorFactory reads from, so this takes effect
# for both loop_interpolation and LoopStructural.interpolators as soon as
# this module (i.e. LoopStructural.modelling.features.fold) has been
# imported once.
interpolator_map[InterpolatorType.DISCRETE_FOLD] = DiscreteFoldInterpolator
