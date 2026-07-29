"""Helpers for fold-normal sign alignment against normal constraints."""

from __future__ import annotations

import numpy as np


_VALID_ALIGNMENT_MODES = {"none", "warn", "correct"}


def resolve_fold_norm_target(
    *,
    fold,
    normal_constraints: np.ndarray,
    fold_norm: float | None,
    dgz_alignment: str | None,
    logger,
    default_fold_norm: float = -1.0,
) -> float:
    """Resolve fold_norm target with optional dgz/normal sign checks.

    Parameters
    ----------
    fold
        Fold event used to evaluate dgz at normal-constraint locations.
    normal_constraints : np.ndarray
        Normal constraints in ``xyz|nx ny nz|w`` format.
    fold_norm : float or None
        Requested fold_norm target. If None, ``default_fold_norm`` is used.
    dgz_alignment : {"none", "warn", "correct"} or None
        Handling mode when dgz appears opposite to normal constraints:
        - ``none``: skip checks.
        - ``warn``: log warning only.
        - ``correct``: log warning and flip target sign.
        ``None`` defaults to ``warn``.
    logger
        Logger used for warnings.
    default_fold_norm : float
        Default fold_norm target when ``fold_norm`` is None.

    Returns
    -------
    float
        Fold normalisation target value to use in ``∇f · dgz = target``.
    """
    target_norm = float(default_fold_norm if fold_norm is None else fold_norm)

    mode = "warn" if dgz_alignment is None else str(dgz_alignment).strip().lower()
    if mode not in _VALID_ALIGNMENT_MODES:
        valid = ", ".join(sorted(_VALID_ALIGNMENT_MODES))
        raise ValueError(f"dgz_alignment must be one of {{{valid}}}, got {dgz_alignment!r}")
    if mode == "none":
        return target_norm

    normals = np.asarray(normal_constraints, dtype=float)
    if normals.ndim != 2 or normals.shape[0] == 0 or normals.shape[1] < 6:
        return target_norm

    points = normals[:, :3]
    normal_vectors = normals[:, 3:6]

    try:
        _, _, dgz = fold.get_deformed_orientation(points)
    except Exception as exc:  # pragma: no cover - defensive fallback
        logger.warning("Could not evaluate dgz for alignment check (%s).", exc)
        return target_norm

    dgz = np.asarray(dgz, dtype=float)
    if dgz.ndim != 2 or dgz.shape[1] != normal_vectors.shape[1]:
        logger.warning(
            "Skipping dgz alignment check due to unexpected dgz shape %s.",
            getattr(dgz, "shape", None),
        )
        return target_norm

    if dgz.shape[0] != normal_vectors.shape[0]:
        n_common = min(dgz.shape[0], normal_vectors.shape[0])
        logger.warning(
            "dgz alignment check received mismatched rows (dgz=%d, normals=%d); using first %d.",
            dgz.shape[0],
            normal_vectors.shape[0],
            n_common,
        )
        if n_common == 0:
            return target_norm
        dgz = dgz[:n_common]
        normal_vectors = normal_vectors[:n_common]

    normal_norm = np.linalg.norm(normal_vectors, axis=1)
    dgz_norm = np.linalg.norm(dgz, axis=1)

    valid = (
        np.all(np.isfinite(normal_vectors), axis=1)
        & np.all(np.isfinite(dgz), axis=1)
        & (normal_norm > 1e-12)
        & (dgz_norm > 1e-12)
    )
    if not np.any(valid):
        return target_norm

    normal_unit = normal_vectors[valid] / normal_norm[valid, None]
    dgz_unit = dgz[valid] / dgz_norm[valid, None]
    median_dot = float(np.median(np.einsum("ij,ij->i", dgz_unit, normal_unit)))

    # A strongly negative median indicates a systematic sign mismatch.
    if median_dot < -0.2:
        logger.warning(
            "Detected dgz opposite to normal constraints (median dot %.3f). "
            "fold_norm may need sign inversion.",
            median_dot,
        )
        if mode == "correct":
            target_norm *= -1.0
            logger.warning(
                "Auto-correcting fold_norm sign to %.3f using dgz_alignment='correct'.",
                target_norm,
            )

    return target_norm
