from LoopStructural.utils.typing import NumericInput
import numpy as np
import numbers
from typing import Tuple


def strikedip2vector(strike: NumericInput, dip: NumericInput) -> np.ndarray:
    """Convert strike and dip to a vector

    Parameters
    ----------
    strike : _type_
        _description_
    dip : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    vec = np.zeros((len(strike), 3))
    s_r = np.deg2rad(strike)
    d_r = np.deg2rad((dip))
    vec[:, 0] = np.sin(d_r) * np.cos(s_r)
    vec[:, 1] = -np.sin(d_r) * np.sin(s_r)
    vec[:, 2] = np.cos(d_r)
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    return vec


def azimuthplunge2vector(
    plunge: NumericInput,
    plunge_dir: NumericInput,
    degrees: bool = True,
) -> np.ndarray:
    """Convert plunge and plunge direction to a vector

    Parameters
    ----------
    plunge : Union[np.ndarray, list]
        array or array like of plunge values
    plunge_dir : Union[np.ndarray, list]
        array or array like of plunge direction values

    Returns
    -------
    np.array
        nx3 vector
    """
    if isinstance(plunge, numbers.Number):
        plunge = np.array([plunge], dtype=float)
    else:
        plunge = np.array(plunge, dtype=float)
    if isinstance(plunge_dir, numbers.Number):
        plunge_dir = np.array([plunge_dir], dtype=float)
    else:
        plunge_dir = np.array(plunge_dir, dtype=float)
    if degrees:
        plunge = np.deg2rad(plunge)
        plunge_dir = np.deg2rad(plunge_dir)
    vec = np.zeros((len(plunge), 3))
    vec[:, 0] = np.sin(plunge_dir) * np.cos(plunge)
    vec[:, 1] = np.cos(plunge_dir) * np.cos(plunge)
    vec[:, 2] = -np.sin(plunge)
    return vec


def normal_vector_to_strike_and_dip(
    normal_vector: NumericInput, degrees: bool = True
) -> np.ndarray:
    """Convert from a normal vector to strike and dip

    Parameters
    ----------
    normal_vector : np.ndarray, list
        array of normal vectors
    degrees : bool, optional
        whether to return in degrees or radians, by default True
    Returns
    -------
    np.ndarray
        2xn array of strike and dip values

    Notes
    ------

    if a 1d array is passed in it is assumed to be a single normal vector
    and cast into a 1x3 array

    """
    normal_vector = np.array(normal_vector)
    if len(normal_vector.shape) == 1:
        normal_vector = normal_vector[None, :]
    # normalise the normal vector
    normal_vector /= np.linalg.norm(normal_vector, axis=1)[:, None]
    dip = np.arccos(normal_vector[:, 2])
    strike = -np.arctan2(normal_vector[:, 1], normal_vector[:, 0])
    if degrees:
        dip = np.rad2deg(dip)
        strike = np.rad2deg(strike)

    return np.array([strike, dip]).T


def rotation(axis: NumericInput, angle: NumericInput) -> np.ndarray:
    """Create a rotation matrix for an axis and angle

    Parameters
    ----------
    axis : Union[np.ndarray, list]
        vector defining the axis of rotation
    angle : Union[np.ndarray, list]
        angle to rotate in degrees

    Returns
    -------
    np.ndarray
        3x3 rotation matrix
    """
    c = np.cos(np.deg2rad(angle))
    s = np.sin((np.deg2rad(angle)))
    C = 1.0 - c
    x = axis[:, 0]
    y = axis[:, 1]
    z = axis[:, 2]
    xs = x * s
    ys = y * s
    zs = z * s
    xC = x * C
    yC = y * C
    zC = z * C
    xyC = x * yC
    yzC = y * zC
    zxC = z * xC
    rotation_mat = np.zeros((axis.shape[0], 3, 3))
    rotation_mat[:, 0, 0] = x * xC + c
    rotation_mat[:, 0, 1] = xyC - zs
    rotation_mat[:, 0, 2] = zxC + ys

    rotation_mat[:, 1, 0] = xyC + zs
    rotation_mat[:, 1, 1] = y * yC + c
    rotation_mat[:, 1, 2] = yzC - xs

    rotation_mat[:, 2, 0] = zxC - ys
    rotation_mat[:, 2, 1] = yzC + xs
    rotation_mat[:, 2, 2] = z * zC + c
    return rotation_mat


def rotate(vector: NumericInput, axis: NumericInput, angle: NumericInput) -> np.ndarray:
    """Rotate a vector about an axis

    Parameters
    ----------
    vector : Union[np.ndarray, list]
        vector to rotate
    alpha : Union[np.ndarray, list]
        axis to rotate about
    beta : Union[np.ndarray, list]
        angle to rotate in degrees

    Returns
    -------
    np.ndarray
        rotated vector
    """
    return np.einsum("ijk,ik->ij", rotation(axis, angle), vector)
    # rotation_mat = rotation(
    #     np.tile(np.array([0, 0, 1])[None, :], (yaw.shape[0], 1)), yaw
    # )
    # vector = np.einsum("ijk,ik->ij", rotation_mat, vector)
    # rotation_mat = rotation(
    #     np.tile(np.array([0, 1, 0])[None, :], (pitch.shape[0], 1)), pitch
    # )
    # vector = np.einsum("ijk,ik->ij", rotation_mat, vector)

    # return vector


def get_vectors(normal: NumericInput) -> Tuple[np.ndarray, np.ndarray]:
    """Find strike and dip vectors for a normal vector.
    Makes assumption the strike vector is horizontal component and the dip is vertical.
    Found by calculating strike and and dip angle and then finding the appropriate vectors

    Parameters
    ----------
    normal : Union[np.ndarray, list]
        input

    Returns
    -------
    np.ndarray, np.ndarray
        strike vector, dip vector
    """
    length = np.linalg.norm(normal, axis=1)[:, None]
    normal /= length  # np.linalg.norm(normal,axis=1)[:,None]
    strikedip = normal_vector_to_strike_and_dip(normal)
    strike_vec = get_strike_vector(strikedip[:, 0])
    strike_vec /= np.linalg.norm(strike_vec, axis=0)[None, :]
    dip_vec = np.cross(strike_vec, normal, axisa=0, axisb=1).T  # (strikedip[:, 0], strikedip[:, 1])
    dip_vec /= np.linalg.norm(dip_vec, axis=0)[None, :]
    return strike_vec * length.T, dip_vec * length.T


def get_strike_vector(strike: NumericInput, degrees: bool = True) -> np.ndarray:
    """Return the vector aligned with the strike direction

    Parameters
    ----------
    strike : np.ndarray
        strike direction
    degrees : bool, optional
        whether to return in degrees or radians, by default True

    Returns
    -------
    np.ndarray
        vector aligned with strike direction

    """
    if isinstance(strike, numbers.Number):
        strike = np.array([strike])
    strike = np.array(strike)
    if degrees:
        strike = np.deg2rad(strike)
    v = np.array(
        [
            np.sin(-strike),
            -np.cos(-strike),
            np.zeros(strike.shape[0]),
        ]
    )

    return v


def get_dip_vector(strike, dip):
    v = np.array(
        [
            -np.cos(np.deg2rad(-strike)) * np.cos(-np.deg2rad(dip)),
            np.sin(np.deg2rad(-strike)) * np.cos(-np.deg2rad(dip)),
            np.sin(-np.deg2rad(dip)),
        ]
    )
    return v


def regular_tetraherdron_for_points(xyz, scale_parameter):
    regular_tetrahedron = np.array(
        [
            [np.sqrt(8 / 9), 0, -1 / 3],
            [-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3],
            [-np.sqrt(2 / 9), -np.sqrt(2 / 3), -1 / 3],
            [0, 0, 1],
        ]
    )
    regular_tetrahedron *= scale_parameter
    tetrahedron = np.zeros((xyz.shape[0], 4, 3))
    tetrahedron[:] = xyz[:, None, :]
    tetrahedron[:, :, :] += regular_tetrahedron[None, :, :]

    return tetrahedron


def gradient_from_tetrahedron(tetrahedron, value):
    """
    Calculate the gradient from a tetrahedron
    """
    tetrahedron = tetrahedron.reshape(-1, 4, 3)
    m = np.array(
        [
            [
                (tetrahedron[:, 1, 0] - tetrahedron[:, 0, 0]),
                (tetrahedron[:, 1, 1] - tetrahedron[:, 0, 1]),
                (tetrahedron[:, 1, 2] - tetrahedron[:, 0, 2]),
            ],
            [
                (tetrahedron[:, 2, 0] - tetrahedron[:, 0, 0]),
                (tetrahedron[:, 2, 1] - tetrahedron[:, 0, 1]),
                (tetrahedron[:, 2, 2] - tetrahedron[:, 0, 2]),
            ],
            [
                (tetrahedron[:, 3, 0] - tetrahedron[:, 0, 0]),
                (tetrahedron[:, 3, 1] - tetrahedron[:, 0, 1]),
                (tetrahedron[:, 3, 2] - tetrahedron[:, 0, 2]),
            ],
        ]
    )
    I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
    m = np.swapaxes(m, 0, 2)
    element_gradients = np.linalg.inv(m)

    element_gradients = element_gradients.swapaxes(1, 2)
    element_gradients = element_gradients @ I
    v = np.sum(element_gradients * value[:, None, :], axis=2)
    return v
