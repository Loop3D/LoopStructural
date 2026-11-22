"""Validation helpers for rigid rotations."""

from __future__ import annotations

import numpy as np


class RigidRotationError(ValueError):
    """Raised when a matrix is not a proper rotation."""


def validate_rigid_rotation(matrix: np.ndarray, dimensions: int, atol: float = 1e-10) -> np.ndarray:
    """Validate and return a copy of a proper rotation matrix.

    Parameters
    ----------
    matrix : np.ndarray
        Candidate rotation matrix.
    dimensions : int
        Required matrix dimensions.
    atol : float, default=1e-10
        Absolute tolerance for orthonormality and determinant checks.

    Returns
    -------
    np.ndarray
        Floating point copy of the validated rotation matrix.

    Raises
    ------
    RigidRotationError
        If the matrix is not a finite proper rotation.
    """
    rotation = np.array(matrix, dtype=float)
    suffix = "; supply a proper rotation (orthonormal, det +1); scale and shear are not supported"

    if rotation.shape != (dimensions, dimensions):
        raise RigidRotationError(f"matrix shape must be ({dimensions}, {dimensions}){suffix}")
    if not np.all(np.isfinite(rotation)):
        raise RigidRotationError(f"matrix entries must be finite{suffix}")

    orthogonality_error = np.max(np.abs(rotation.T @ rotation - np.eye(dimensions)))
    if orthogonality_error > atol:
        raise RigidRotationError(f"matrix must be orthonormal{suffix}")

    determinant_error = abs(np.linalg.det(rotation) - 1.0)
    if determinant_error > atol:
        raise RigidRotationError(f"matrix determinant must be +1{suffix}")

    return rotation
