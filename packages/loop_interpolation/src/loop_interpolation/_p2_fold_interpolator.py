"""
Piecewise quadratic interpolator using fold constraints.

Summary of findings_1_1.md (what exists / what was added)
---------------------------------------------------------
Existing:
- ``P2Interpolator`` and the P2 supports. ``evaluate_shape_derivatives``
  gives dN of shape (n, 3, 10), node order [v0..v3, e23, e03, e01, e12, e13, e02].
- ``setup_with_fold_constraints`` (_fold_setup.py) and
  ``resolve_fold_norm_target`` (_fold_norm_alignment.py), both reused as is.
- ``add_constraints_to_least_squares`` drops NaN rows, normalises the rows of A,
  then applies the weight ``w``. A per-row ``w`` is therefore applied after
  normalisation.
- ``get_quadrature_points`` on the support is a triangle (face) rule for edge
  jumps. It does not fit a volume rule.
Added (Phase 1):
- ``_quadrature.py``: ``tet_quadrature`` / ``tet_quadrature_points``
  (centroid, gauss4, weights sum to 1, tet index kept).
- P2 ``evaluate_shape_derivatives`` returns (dN, elements, inside) with NaN
  rows for outside points.
- This class. ``copy()`` is overridden because the base returns P2Interpolator.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from loop_common.logging import get_logger as getLogger

from ._fold_event import FoldEvent
from ._fold_norm_alignment import resolve_fold_norm_target
from ._fold_setup import setup_with_fold_constraints
from ._interpolatortype import InterpolatorType
from ._p2interpolator import P2Interpolator
from ._quadrature import tet_quadrature_points

logger = getLogger(__name__)

_SAMPLE_RULES = ("centroid", "gauss4")


class P2FoldInterpolator(P2Interpolator):
    """P2 interpolator with fold constraints (Laurent et al., 2016).

    ``nelements`` (when built through ``InterpolatorFactory``) is the number of
    tetrahedra. P2 has a node on every edge, so 1e4 gives about 56k tets and 86k
    unknowns, and 1e5 gives about 526k tets and 765k unknowns (solve about 2.5 min).
    The mesh is larger than the model box because of the buffer.
    """

    def __init__(self, support, fold: FoldEvent | None = None):
        P2Interpolator.__init__(self, support)
        self.type = InterpolatorType.PIECEWISE_QUADRATIC_FOLD
        self.fold = fold
        logger.info(
            f"P2Fold support: {support.elements.shape[0]} tets, {support.nodes.shape[0]} unknowns. "
            "nelements counts tets for P2, so cost grows quickly with it"
        )

    def copy(self):
        return P2FoldInterpolator(self.support, fold=self.fold)

    def setup_interpolator(self, **kwargs):
        report = setup_with_fold_constraints(
            fold=self.fold,
            kwargs=kwargs,
            interpolator_name=self.__class__.__name__,
            base_setup=super().setup_interpolator,
            add_fold_constraints=self.add_fold_constraints,
            finalize_report=self.finalize_setup_diagnostics_report,
        )
        # set after the call: includes fold_weights, which the base setup does not see
        self._setup_kwargs = dict(kwargs)
        return report

    def add_fold_constraints(
        self,
        fold_orientation=10.0,
        fold_axis_w=10.0,
        fold_normalisation=1.0,
        fold_norm=-1.0,
        dgz_alignment="warn",
        sample_rule="gauss4",
        mask_fn: Callable | None = None,
    ):
        """Add the three fold constraint families at quadrature points.

        Parameters
        ----------
        fold_orientation, fold_axis_w, fold_normalisation : float or None
            weights of the families. None skips a family.
        fold_norm : float
            norm target, see ``resolve_fold_norm_target``
        dgz_alignment : {"none", "warn", "correct"}
        sample_rule : {"gauss4", "centroid"}
            quadrature rule used to place sample points in every tetrahedron
        mask_fn : callable, optional
            called with the sample points; points where it returns True get zero weight
        """
        if sample_rule not in _SAMPLE_RULES:
            raise ValueError(f"Unknown sample_rule '{sample_rule}', use one of {_SAMPLE_RULES}")
        support = self.support
        points, tet_idx, q = tet_quadrature_points(support.nodes, support.elements, sample_rule)
        dN, _out_elements, inside = support.evaluate_shape_derivatives(points, elements=tet_idx)
        deformed_orientation, fold_axis, dgz = self.fold.get_deformed_orientation(points)
        # energy per tet must not depend on the rule: weight the residual by w*sqrt(q)
        sq = np.sqrt(q)
        if mask_fn is not None:
            sq = sq * (~np.asarray(mask_fn(points), dtype=bool))
        # keep finite rows only, so A, B, idc and w stay aligned
        valid = np.asarray(inside, dtype=bool) & np.all(np.isfinite(dN), axis=(1, 2))
        idc_all = support.elements[tet_idx]

        def _add(vectors, target, w, name):
            v = valid & np.all(np.isfinite(vectors), axis=1)
            # no volume factor: row normalisation would cancel it
            A = np.einsum("ij,ijk->ik", vectors[v], dN[v])
            B = np.full(A.shape[0], float(target))
            self.add_constraints_to_least_squares(A, B, idc_all[v], w=float(w) * sq[v], name=name)

        if fold_orientation is not None:
            logger.info(f"Adding fold orientation constraint to w = {fold_orientation}")
            _add(deformed_orientation, 0.0, fold_orientation, "fold orientation")
        if fold_axis_w is not None:
            logger.info(f"Adding fold axis constraint to w = {fold_axis_w}")
            _add(fold_axis, 0.0, fold_axis_w, "fold axis")
        if fold_normalisation is not None:
            logger.info(f"Adding fold normalisation constraint to w = {fold_normalisation}")
            target_norm = resolve_fold_norm_target(
                fold=self.fold,
                normal_constraints=self.get_norm_constraints(),
                fold_norm=fold_norm,
                dgz_alignment=dgz_alignment,
                logger=logger,
                default_fold_norm=-1.0,
            )
            # Decision: unlike DiscreteFoldInterpolator, B is NOT multiplied by
            # fold_normalisation (nor by volume). That scaling changes the target
            # norm with the weight, which is a quirk. Here the weight only sets
            # the importance and the target stays fold_norm.
            _add(dgz, target_norm, fold_normalisation, "fold normalisation")
