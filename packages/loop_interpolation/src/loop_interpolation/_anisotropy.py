"""Locally varying (anisotropic) regularisation constraints.

These are plain functions, not interpolator subclasses: they add a bundle of
constraints (two gradient-orthogonality constraints, a normal-magnitude
constraint, and directional regularisation) to an *existing*
:class:`FiniteDifferenceInterpolator` or :class:`P1Interpolator` instance,
built entirely out of the generic constraint machinery those classes already
expose (``add_gradient_orthogonal_constraints``, ``add_directional_regularisation``).
The only thing genuinely specific to "anisotropy" here is evaluating a
direction provider (an object with ``get_directions(points) ->
(primary, secondary, normal)``) at the right points and threading the results
through those generic calls.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from loop_common.logging import get_logger as getLogger
from loop_common.math import rng

from ._regularisation import DirectionalRegularisation

logger = getLogger(__name__)

_VALID_ALIGNMENT_MODES = {"none", "warn", "correct"}
_DEFAULT_REGULARISATION_WEIGHTS = (0.1, 0.01, 0.01)


def resolve_anisotropy_norm_target(
    *,
    anisotropy,
    normal_constraints: np.ndarray,
    norm_target: float | None,
    norm_alignment: str | None,
    logger,
    default_norm_target: float = -1.0,
) -> float:
    """Resolve the normal-direction norm target with optional sign checks.

    Parameters
    ----------
    anisotropy
        Anisotropy direction provider used to evaluate the normal direction
        at normal-constraint locations.
    normal_constraints : np.ndarray
        Normal constraints in ``xyz|nx ny nz|w`` format.
    norm_target : float or None
        Requested norm target. If None, ``default_norm_target`` is used.
    norm_alignment : {"none", "warn", "correct"} or None
        Handling mode when the anisotropy normal direction appears opposite
        to normal constraints:
        - ``none``: skip checks.
        - ``warn``: log warning only.
        - ``correct``: log warning and flip target sign.
        ``None`` defaults to ``warn``.
    logger
        Logger used for warnings.
    default_norm_target : float
        Default norm target when ``norm_target`` is None.

    Returns
    -------
    float
        Norm target value to use in ``∇f · normal = target``.
    """
    target_norm = float(default_norm_target if norm_target is None else norm_target)

    mode = "warn" if norm_alignment is None else str(norm_alignment).strip().lower()
    if mode not in _VALID_ALIGNMENT_MODES:
        valid = ", ".join(sorted(_VALID_ALIGNMENT_MODES))
        raise ValueError(f"norm_alignment must be one of {{{valid}}}, got {norm_alignment!r}")
    if mode == "none":
        return target_norm

    normals = np.asarray(normal_constraints, dtype=float)
    if normals.ndim != 2 or normals.shape[0] == 0 or normals.shape[1] < 6:
        return target_norm

    points = normals[:, :3]
    normal_vectors = normals[:, 3:6]

    try:
        _, _, anisotropy_normal = anisotropy.get_directions(points)
    except (AttributeError, TypeError, ValueError, RuntimeError) as exc:  # pragma: no cover - defensive fallback
        logger.warning("Could not evaluate anisotropy normal for alignment check (%s).", exc)
        return target_norm

    anisotropy_normal = np.asarray(anisotropy_normal, dtype=float)
    if anisotropy_normal.ndim != 2 or anisotropy_normal.shape[1] != normal_vectors.shape[1]:
        logger.warning(
            "Skipping normal alignment check due to unexpected shape %s.",
            getattr(anisotropy_normal, "shape", None),
        )
        return target_norm

    if anisotropy_normal.shape[0] != normal_vectors.shape[0]:
        n_common = min(anisotropy_normal.shape[0], normal_vectors.shape[0])
        logger.warning(
            "Normal alignment check received mismatched rows (anisotropy=%d, normals=%d); using first %d.",
            anisotropy_normal.shape[0],
            normal_vectors.shape[0],
            n_common,
        )
        if n_common == 0:
            return target_norm
        anisotropy_normal = anisotropy_normal[:n_common]
        normal_vectors = normal_vectors[:n_common]

    normal_norm = np.linalg.norm(normal_vectors, axis=1)
    anisotropy_norm = np.linalg.norm(anisotropy_normal, axis=1)

    valid = (
        np.all(np.isfinite(normal_vectors), axis=1)
        & np.all(np.isfinite(anisotropy_normal), axis=1)
        & (normal_norm > 1e-12)
        & (anisotropy_norm > 1e-12)
    )
    if not np.any(valid):
        return target_norm

    normal_unit = normal_vectors[valid] / normal_norm[valid, None]
    anisotropy_unit = anisotropy_normal[valid] / anisotropy_norm[valid, None]
    median_dot = float(np.median(np.einsum("ij,ij->i", anisotropy_unit, normal_unit)))

    # A strongly negative median indicates a systematic sign mismatch.
    if median_dot < -0.2:
        logger.warning(
            "Detected anisotropy normal opposite to normal constraints (median dot %.3f). "
            "norm_target may need sign inversion.",
            median_dot,
        )
        if mode == "correct":
            target_norm *= -1.0
            logger.warning(
                "Auto-correcting norm_target sign to %.3f using norm_alignment='correct'.",
                target_norm,
            )

    return target_norm


def add_fd_anisotropy_constraints(
    interpolator,
    anisotropy,
    primary_w: float | None = 10.0,
    secondary_w: float | None = 10.0,
    regularisation_weights=_DEFAULT_REGULARISATION_WEIGHTS,
    normal_w: float | None = 1.0,
    norm_target: float | None = -1.0,
    norm_alignment: str = "warn",
    mask_fn: Callable | None = None,
) -> None:
    """Add locally varying anisotropy constraints to a finite difference interpolator.

    Direction vectors are evaluated at every grid *node* by calling
    ``anisotropy.get_directions(interpolator.support.nodes)``. The returned
    vectors are then used to build, entirely via the interpolator's existing
    generic constraint methods:

    * Gradient-orientation constraints (dot product = 0) for the primary
      and secondary directions.
    * A norm constraint forcing the gradient magnitude in the normal
      direction to be ``norm_target``.
    * Anisotropic regularisation using directional second derivatives,
      applied with different weights for the three directions.

    Parameters
    ----------
    interpolator
        A :class:`FiniteDifferenceInterpolator` (or compatible) instance to
        add constraints to.
    anisotropy
        Direction provider with a ``get_directions(points)`` method
        returning ``(primary, secondary, normal)`` unit-vector fields (each
        an ``(n_points, 3)`` array).
    primary_w : float or None
        Weight for the constraint that :math:`\\nabla f` is perpendicular
        to the primary direction. Pass ``None`` to skip.
    secondary_w : float or None
        Weight for the constraint that :math:`\\nabla f` is perpendicular
        to the secondary direction. Pass ``None`` to skip.
    regularisation_weights : list of three floats or None
        Weights ``[w0, w1, w2]`` for the anisotropic regularisation applied
        in the normal, primary, and secondary directions respectively.
        Defaults to ``(0.1, 0.01, 0.01)``. Pass ``None`` to skip entirely.
    normal_w : float or None
        Weight for the norm constraint. Pass ``None`` to skip.
    norm_target : float or None
        Target gradient projection in the normal direction. Defaults to
        ``-1.0``.
    norm_alignment : {"none", "warn", "correct"}
        How to handle detected sign mismatch between the normal direction
        and normal constraints (if present).
    mask_fn : callable or None
        If provided, called with ``(n_nodes, 3)`` node positions; nodes for
        which the function returns ``True`` are excluded from all
        anisotropy constraints.
    """
    if anisotropy is None:
        raise ValueError("anisotropy must not be None")

    node_positions = interpolator.support.nodes  # (n_nodes, 3)
    primary, secondary, normal = anisotropy.get_directions(node_positions)

    weight = np.ones(interpolator.support.n_nodes, dtype=float)
    if mask_fn is not None:
        weight[mask_fn(node_positions)] = 0.0
    active = weight > 0.0
    active_pos = node_positions[active]

    if primary_w is not None:
        logger.info(f"Adding primary orientation constraint  w = {primary_w}")
        interpolator.add_gradient_orthogonal_constraints(
            active_pos, primary[active], w=primary_w, name="anisotropy primary"
        )

    if secondary_w is not None:
        logger.info(f"Adding secondary orientation constraint  w = {secondary_w}")
        interpolator.add_gradient_orthogonal_constraints(
            active_pos, secondary[active], w=secondary_w, name="anisotropy secondary"
        )

    if normal_w is not None:
        target_norm = resolve_anisotropy_norm_target(
            anisotropy=anisotropy,
            normal_constraints=interpolator.get_norm_constraints(),
            norm_target=norm_target,
            norm_alignment=norm_alignment,
            logger=logger,
            default_norm_target=-1.0,
        )
        logger.info(f"Adding normal magnitude constraint  w = {normal_w}")
        interpolator.add_gradient_orthogonal_constraints(
            active_pos, normal[active], w=normal_w, b=target_norm, name="anisotropy normal"
        )

    if regularisation_weights is not None:
        logger.info(
            f"Adding anisotropic regularisation  w = {regularisation_weights[0]}, "
            f"{regularisation_weights[1]}, {regularisation_weights[2]}"
        )

        def _masked_direction(vectors):
            masked = np.asarray(vectors, dtype=float).copy()
            masked[~active] = 0.0
            return masked

        interpolator.add_directional_regularisation(
            (
                DirectionalRegularisation(
                    weight=regularisation_weights[0],
                    direction=_masked_direction(normal),
                    name="anisotropy regularisation 1",
                ),
                DirectionalRegularisation(
                    weight=regularisation_weights[1],
                    direction=_masked_direction(primary),
                    name="anisotropy regularisation 2",
                ),
                DirectionalRegularisation(
                    weight=regularisation_weights[2],
                    direction=_masked_direction(secondary),
                    name="anisotropy regularisation 3",
                ),
            )
        )


def add_element_anisotropy_constraints(
    interpolator,
    anisotropy,
    primary_w: float | None = 10.0,
    secondary_w: float | None = 10.0,
    regularisation_weights=_DEFAULT_REGULARISATION_WEIGHTS,
    normal_w: float | None = 1.0,
    norm_target: float | None = -1.0,
    norm_alignment: str = "warn",
    step: int = 2,
    mask_fn: Callable | None = None,
) -> None:
    """Add locally varying anisotropy constraints to an element-based (mesh) interpolator.

    Works with any interpolator built on a mesh support exposing
    ``support.barycentre``/``support.n_elements`` and the shared
    ``add_gradient_orthogonal_constraints``/``add_directional_regularisation``
    API — currently :class:`P1Interpolator` and :class:`P2Interpolator`.
    (:class:`FiniteDifferenceInterpolator`'s node-based grid uses
    :func:`add_fd_anisotropy_constraints` instead.)

    Direction vectors are evaluated at every mesh element's barycentre by
    calling ``anisotropy.get_directions(interpolator.support.barycentre)``.
    A random ``step``-subsample of elements is used for each of the three
    gradient-orientation/norm constraint families (independently reshuffled
    per family, matching the historical fold-interpolator behaviour), while
    the directional regularisation term is applied densely (masked but not
    subsampled), evaluated lazily wherever the underlying edge-jump
    regularisation needs it.

    Parameters
    ----------
    interpolator
        A :class:`P1Interpolator`, :class:`P2Interpolator` (or compatible)
        instance to add constraints to.
    anisotropy
        Direction provider with a ``get_directions(points)`` method
        returning ``(primary, secondary, normal)`` unit-vector fields.
    primary_w : double
        weight for the primary direction constraint in the least squares system
    secondary_w : double
        weight for the secondary direction constraint in the least squares system
    regularisation_weights : list
        weight for the anisotropic regularisation in the least squares system
    normal_w : double
        weight for the normal-direction norm constraint in the least squares system
    norm_target
        length of the interpolation norm in the least squares system
    norm_alignment : {"none", "warn", "correct"}
        How to handle detected sign mismatch between the normal direction
        and normal constraints (if present).
    step : int
        array step for adding constraints
    mask_fn : callable or None
        If provided, called with ``(n_elements, 3)`` barycentre positions;
        elements for which the function returns ``True`` are excluded.

    Notes
    -----
    For more information about the underlying constraints see EPSL paper by Gautier Laurent 2016
    """
    if anisotropy is None:
        raise ValueError("anisotropy must not be None")

    support = interpolator.support
    barycentre = support.barycentre
    primary, secondary, normal = anisotropy.get_directions(barycentre)

    weight = np.ones(support.n_elements, dtype=float)
    if mask_fn is not None:
        weight[mask_fn(barycentre)] = 0.0

    element_idx = np.arange(support.n_elements)
    rng.shuffle(element_idx)

    def _active_subsample():
        rng.shuffle(element_idx)
        selected = element_idx[::step]
        return selected[weight[selected] > 0.0]

    if primary_w is not None:
        logger.info(f"Adding primary orientation constraint  w = {primary_w}")
        active = _active_subsample()
        interpolator.add_gradient_orthogonal_constraints(
            barycentre[active], primary[active], w=primary_w, name="anisotropy primary"
        )

    if secondary_w is not None:
        logger.info(f"Adding secondary orientation constraint  w = {secondary_w}")
        active = _active_subsample()
        interpolator.add_gradient_orthogonal_constraints(
            barycentre[active], secondary[active], w=secondary_w, name="anisotropy secondary"
        )

    if normal_w is not None:
        target_norm = resolve_anisotropy_norm_target(
            anisotropy=anisotropy,
            normal_constraints=interpolator.get_norm_constraints(),
            norm_target=norm_target,
            norm_alignment=norm_alignment,
            logger=logger,
            default_norm_target=-1.0,
        )
        logger.info(f"Adding normal magnitude constraint  w = {normal_w}")
        active = _active_subsample()
        interpolator.add_gradient_orthogonal_constraints(
            barycentre[active], normal[active], w=normal_w, b=target_norm, name="anisotropy normal"
        )

    if regularisation_weights is not None:
        logger.info(
            f"Adding anisotropic regularisation  w = {regularisation_weights[0]}, "
            f"{regularisation_weights[1]}, {regularisation_weights[2]}"
        )

        def _masked_direction(component_index: int):
            def _provider(points: np.ndarray) -> np.ndarray:
                primary_dir, secondary_dir, normal_dir = anisotropy.get_directions(points)
                vectors = (normal_dir, primary_dir, secondary_dir)[component_index]
                if mask_fn is not None:
                    masked = np.asarray(vectors, dtype=float).copy()
                    masked[mask_fn(points)] = 0.0
                    return masked
                return vectors

            return _provider

        interpolator.add_directional_regularisation(
            (
                DirectionalRegularisation(
                    weight=regularisation_weights[0],
                    direction=_masked_direction(0),
                    name="anisotropy regularisation 1",
                ),
                DirectionalRegularisation(
                    weight=regularisation_weights[1],
                    direction=_masked_direction(1),
                    name="anisotropy regularisation 2",
                ),
                DirectionalRegularisation(
                    weight=regularisation_weights[2],
                    direction=_masked_direction(2),
                    name="anisotropy regularisation 3",
                ),
            )
        )
