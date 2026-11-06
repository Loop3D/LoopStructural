"""Convert between vectors and orientation angles.

Convention: x is east, y is north, z is up. All angles are in degrees.
A plane has an upward normal. A lineation is an axis that points down.
"""

from __future__ import annotations

import numpy as np

from loop_link.errors import InvalidInputError

__all__ = [
    "directions_to_lineation",
    "lineation_to_directions",
    "normals_to_orientation",
    "orientation_to_normals",
]


def _unit_vectors(vectors: object, name: str) -> np.ndarray:
    """Return the input as unit vectors with shape (n, 3).

    Parameters
    ----------
    vectors : array-like
        Vectors with shape (n, 3) or (3,).
    name : str
        The name of the input. The name is used in error messages.

    Returns
    -------
    numpy.ndarray
        Unit vectors with shape (n, 3).

    Raises
    ------
    InvalidInputError
        If the shape is wrong, or a vector is zero-length or not finite.
    """
    try:
        arr = np.asarray(vectors, dtype=float)
    except (TypeError, ValueError) as exc:
        raise InvalidInputError(f"The {name} must be numbers. Check the input.") from exc
    arr = np.atleast_2d(arr)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise InvalidInputError(f"The {name} must have shape (n, 3) or (3,). Check the input.")
    if not np.all(np.isfinite(arr)):
        raise InvalidInputError(f"The {name} must be finite numbers. Check the input.")
    lengths = np.linalg.norm(arr, axis=1)
    if not np.all(np.isfinite(lengths)) or np.any(lengths == 0.0):
        raise InvalidInputError(f"The {name} must not have a zero-length vector. Check the input.")
    return arr / lengths[:, None]


def _angles(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the azimuth of (x, y) in degrees, in [0, 360).

    The azimuth is 0.0 when x and y are both zero.
    """
    flat = (x == 0.0) & (y == 0.0)
    # Adding 0.0 turns -0.0 into 0.0.
    azimuth = np.degrees(np.arctan2(np.where(flat, 0.0, x), np.where(flat, 1.0, y))) % 360.0 + 0.0
    # A tiny negative angle can round to 360.0 after the modulo.
    azimuth = np.where(azimuth >= 360.0, 0.0, azimuth)
    return azimuth


def normals_to_orientation(normals: object) -> tuple[np.ndarray, np.ndarray]:
    """Convert plane normals to dip direction and dip.

    A downward normal is flipped to point up before the conversion.

    Parameters
    ----------
    normals : array-like
        Normal vectors with shape (n, 3) or (3,).

    Returns
    -------
    dip_direction : numpy.ndarray
        Dip direction in degrees, in [0, 360).
    dip : numpy.ndarray
        Dip in degrees, in [0, 90].

    Raises
    ------
    InvalidInputError
        If a vector is zero-length or not finite.
    """
    n = _unit_vectors(normals, "normals")
    n = np.where(n[:, 2:3] < 0.0, -n, n)
    dip = np.degrees(np.arccos(np.clip(n[:, 2], -1.0, 1.0)))
    dip_direction = _angles(n[:, 0], n[:, 1])
    return dip_direction, dip


def directions_to_lineation(directions: object) -> tuple[np.ndarray, np.ndarray]:
    """Convert line directions to trend and plunge.

    An upward direction is flipped to point down before the conversion.

    Parameters
    ----------
    directions : array-like
        Direction vectors with shape (n, 3) or (3,).

    Returns
    -------
    trend : numpy.ndarray
        Trend in degrees, in [0, 360).
    plunge : numpy.ndarray
        Plunge in degrees, in [0, 90].

    Raises
    ------
    InvalidInputError
        If a vector is zero-length or not finite.
    """
    d = _unit_vectors(directions, "directions")
    d = np.where(d[:, 2:3] > 0.0, -d, d)
    plunge = np.degrees(np.arcsin(np.clip(-d[:, 2], -1.0, 1.0)))
    trend = _angles(d[:, 0], d[:, 1])
    return trend, plunge


def orientation_to_normals(dip_direction: object, dip: object) -> np.ndarray:
    """Convert dip direction and dip to upward unit normals.

    Parameters
    ----------
    dip_direction : array-like
        Dip direction in degrees.
    dip : array-like
        Dip in degrees.

    Returns
    -------
    numpy.ndarray
        Unit normals with shape (n, 3).
    """
    a = np.radians(np.atleast_1d(np.asarray(dip_direction, dtype=float)))
    d = np.radians(np.atleast_1d(np.asarray(dip, dtype=float)))
    return np.column_stack((np.sin(d) * np.sin(a), np.sin(d) * np.cos(a), np.cos(d) + 0.0 * a))


def lineation_to_directions(trend: object, plunge: object) -> np.ndarray:
    """Convert trend and plunge to downward unit direction vectors.

    Parameters
    ----------
    trend : array-like
        Trend in degrees.
    plunge : array-like
        Plunge in degrees.

    Returns
    -------
    numpy.ndarray
        Unit direction vectors with shape (n, 3).
    """
    t = np.radians(np.atleast_1d(np.asarray(trend, dtype=float)))
    p = np.radians(np.atleast_1d(np.asarray(plunge, dtype=float)))
    return np.column_stack((np.cos(p) * np.sin(t), np.cos(p) * np.cos(t), -np.sin(p) + 0.0 * t))
