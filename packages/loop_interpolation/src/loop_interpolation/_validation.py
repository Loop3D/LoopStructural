"""Input validation helpers for loop_interpolation.

This module provides centralized validation for constraint inputs to ensure
consistent error handling, clear error messages, and early detection of
invalid constraint combinations.

NaN handling contract
---------------------
- **Weight column NaN** → silently replaced with 1.0 (default weight).
  This lets callers pass ``np.nan`` as a sentinel for "use default weight",
  which is a common pattern when assembling constraints from DataFrames.
- **Position / data column NaN or inf** → rows are silently dropped with a
  warning.  This preserves backward-compatible behaviour where constraints
  that fall outside the model or contain missing values are ignored.
"""

from __future__ import annotations

import logging

import numpy as np

_logger = logging.getLogger(__name__)


def _drop_nan_data_rows(points: np.ndarray, data_cols: slice, name: str) -> np.ndarray:
    """Drop rows that contain NaN or inf in the position / data columns.

    Weight columns are intentionally excluded from this check because NaN
    weights are handled separately.

    Parameters
    ----------
    points : np.ndarray
        Constraint array, already converted to float64.
    data_cols : slice
        Slice selecting the position / data columns (excluding weight).
    name : str
        Human-readable constraint type name used in the warning message.

    Returns
    -------
    np.ndarray
        Array with invalid rows removed; the original array is returned
        unchanged when no rows are dropped.
    """
    bad_rows = ~np.isfinite(points[:, data_cols]).all(axis=1)
    if np.any(bad_rows):
        n_dropped = int(np.sum(bad_rows))
        _logger.warning(
            "%d %s row(s) dropped: position or data columns contain NaN or inf. "
            "Provide finite values to include them.",
            n_dropped,
            name,
        )
        points = points[~bad_rows]
    return points


class ValidationError(ValueError):
    """Base exception for constraint validation errors."""


class ShapeError(ValidationError):
    """Exception for shape mismatches in constraint arrays."""


class DtypeError(ValidationError):
    """Exception for data type mismatches in constraint arrays."""


class FiniteValueError(ValidationError):
    """Exception for non-finite values in constraint arrays."""


class VectorError(ValidationError):
    """Exception for invalid vector/direction constraints."""


class WeightError(ValidationError):
    """Exception for invalid weight values."""


class UnsupportedCombinationError(ValidationError):
    """Exception for unsupported constraint combinations."""


def _ensure_float_array(arr: np.ndarray, name: str = "array") -> np.ndarray:
    """Convert array to float type with error handling.

    Parameters
    ----------
    arr : np.ndarray
        Input array to convert
    name : str
        Name of array for error messages

    Returns
    -------
    np.ndarray
        Array converted to float64

    Raises
    ------
    DtypeError
        If array cannot be converted to float
    """
    try:
        return np.asarray(arr, dtype=np.float64)
    except (TypeError, ValueError) as e:
        raise DtypeError(
            f"{name} could not be converted to float64. "
            f"All constraint arrays must contain numeric values. Error: {e}"
        )


def _check_finite(arr: np.ndarray, name: str = "array") -> None:
    """Validate all values in array are finite (not NaN or inf).

    Parameters
    ----------
    arr : np.ndarray
        Array to validate
    name : str
        Name of array for error messages

    Raises
    ------
    FiniteValueError
        If any non-finite values are found
    """
    non_finite = ~np.isfinite(arr)
    if np.any(non_finite):
        n_invalid = int(np.sum(non_finite))
        invalid_indices = np.where(non_finite)
        raise FiniteValueError(
            f"{name} contains {n_invalid} non-finite values (NaN or inf). "
            f"All constraint values must be finite. "
            f"Found at indices: {invalid_indices}"
        )


def check_unsupported_combinations(
    data: dict,
    constraint_type: str,
) -> None:
    """Check for unsupported constraint combinations.

    Parameters
    ----------
    data : dict
        Dictionary of all constraints
    constraint_type : str
        Type of constraint being added

    Raises
    ------
    UnsupportedCombinationError
        If unsupported combinations are detected
    """
    supported_constraint_types = {
        "value",
        "gradient",
        "normal",
        "tangent",
        "interface",
        "inequality",
        "inequality_pairs",
    }

    if constraint_type not in supported_constraint_types:
        raise UnsupportedCombinationError(
            f"Unsupported constraint type '{constraint_type}'. "
            f"Supported types are: {sorted(supported_constraint_types)}"
        )

    # Future extension point for solver-specific constraint incompatibilities.
    if not isinstance(data, dict):
        raise UnsupportedCombinationError(
            "Constraint container must be a dict mapping constraint families to arrays."
        )
