"""Shared setup helpers for fold-aware interpolators."""

from __future__ import annotations

import numpy as np


def setup_with_fold_constraints(
    *,
    fold,
    kwargs: dict,
    interpolator_name: str,
    base_setup,
    add_fold_constraints,
    finalize_report,
):
    """Run common fold interpolator setup sequence.

    Parameters
    ----------
    fold
        Fold event instance or None.
    kwargs : dict
        Setup kwargs passed to setup_interpolator.
    interpolator_name : str
        Name used in error messages.
    base_setup : callable
        Parent setup_interpolator callable.
    add_fold_constraints : callable
        Method that applies fold-specific constraints.
    finalize_report : callable
        Callable that returns final diagnostics report.
    """
    if fold is None:
        raise RuntimeError(
            f"{interpolator_name}: no fold event set. Assign self.fold before calling setup_interpolator."
        )

    setup_kwargs = dict(kwargs)
    fold_weights = setup_kwargs.pop("fold_weights", {})
    base_setup(**setup_kwargs)
    add_fold_constraints(**fold_weights)
    return finalize_report()


def sample_fold_geometry(
    interpolator, fold, local_points: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample fold geometry at local model points.

    Parameters
    ----------
    interpolator
        Interpolator that maps local points to world coordinates.
    fold
        Fold event that provides ``get_deformed_orientation``. The fold
        provider takes world points and returns world unit directions.
    local_points : np.ndarray
        Points in local model coordinates, shape (n, 3).

    Returns
    -------
    tuple of np.ndarray
        ``(world_points, deformed, axis, dgz)``. Each array is float.
        NaN rows in the fold output stay NaN.
    """
    world_points = interpolator.local_to_world_points(np.asarray(local_points, dtype=float))
    deformed, axis, dgz = fold.get_deformed_orientation(world_points)
    return (
        world_points,
        np.asarray(deformed, dtype=float),
        np.asarray(axis, dtype=float),
        np.asarray(dgz, dtype=float),
    )


def world_norm_constraints(interpolator) -> np.ndarray:
    """Return norm constraints with points and vectors in world coordinates.

    Parameters
    ----------
    interpolator
        Interpolator that provides ``get_norm_constraints`` in local
        coordinates and maps local points and vectors to world coordinates.

    Returns
    -------
    np.ndarray
        Constraint array. Columns 0:3 are world points and columns 3:6
        are world vectors. Arrays with fewer than six columns, no rows,
        or other than two dimensions are returned unchanged.
    """
    n = np.asarray(interpolator.get_norm_constraints(), dtype=float)
    if n.ndim != 2 or n.shape[0] == 0 or n.shape[1] < 6:
        return n
    out = n.copy()
    out[:, 0:3] = interpolator.local_to_world_points(n[:, 0:3])
    out[:, 3:6] = interpolator.local_to_world_vectors(n[:, 3:6])
    return out
