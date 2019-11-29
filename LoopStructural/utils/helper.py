import logging

import numpy as np

logger = logging.getLogger(__name__)


def plunge_and_plunge_dir_to_vector(plunge, plunge_dir):
    plunge = np.deg2rad(plunge)
    plunge_dir = np.deg2rad(plunge_dir)
    vec = np.zeros(3)
    vec[0] = np.sin(plunge_dir) * np.cos(plunge)
    vec[1] = np.cos(plunge_dir) * np.cos(plunge)
    vec[2] = -np.sin(plunge)
    return vec


def create_surface(bounding_box, nstep):
    x = np.linspace(bounding_box[0, 0], bounding_box[1, 0], nstep[0])  #
    y = np.linspace(bounding_box[0, 1], bounding_box[1, 1], nstep[1])
    xx, yy = np.meshgrid(x, y, indexing='xy')

    def gi(i, j):
        return i + j * nstep[0]

    corners = np.array([[0, 1, 0, 1], [0, 0, 1, 1]])
    i = np.arange(0, nstep[0] - 1)

    j = np.arange(0, nstep[1] - 1)
    ii, jj = np.meshgrid(i, j, indexing='ij')
    corner_gi = gi(ii[:, :, None] + corners[None, None, 0, :, ],
                   jj[:, :, None] + corners[None, None, 1, :, ])
    corner_gi = corner_gi.reshape((nstep[0] - 1) * (nstep[1] - 1), 4)
    tri = np.vstack([corner_gi[:, :3], corner_gi[:, 1:]])
    return tri, xx.flatten(), yy.flatten()


def get_vectors(normal):
    strikedip = normal_vector_to_strike_and_dip(normal)
    dip_vec = get_dip_vector(strikedip[:, 0], strikedip[:, 1])
    strike_vec = get_strike_vector(strikedip[:, 0])
    return strike_vec, dip_vec


def get_strike_vector(strike):
    """

    Parameters
    ----------
    strike

    Returns
    -------

    """
    v = np.array([-np.sin(np.deg2rad(-strike)),
                  np.cos(np.deg2rad(-strike)),
                  np.zeros(strike.shape[0])
                  ])
    return v


def get_dip_vector(strike, dip):
    v = np.array([-np.cos(np.deg2rad(-strike)) * np.cos(-np.deg2rad(dip)),
                  np.sin(np.deg2rad(-strike)) * np.cos(-np.deg2rad(dip)),
                  np.sin(-np.deg2rad(dip))
                  ])
    return v


def rotation(axis, angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin((np.deg2rad(angle)))
    C = 1.0 - c
    x = axis[0]
    y = axis[1]
    z = axis[2]
    xs = x * s
    ys = y * s
    zs = z * s
    xC = x * C
    yC = y * C
    zC = z * C
    xyC = x * yC
    yzC = y * zC
    zxC = z * xC
    rotation_mat = np.zeros((3, 3))
    rotation_mat[0][0] = x * xC + c
    rotation_mat[0][1] = xyC - zs
    rotation_mat[0][2] = zxC + ys

    rotation_mat[1][0] = xyC + zs
    rotation_mat[1][1] = y * yC + c
    rotation_mat[1][2] = yzC - xs

    rotation_mat[2][0] = zxC - ys
    rotation_mat[2][1] = yzC + xs
    rotation_mat[2][2] = z * zC + c
    return rotation_mat


def normalz(gx):
    gxn = (2. / (np.max(gx[~np.isnan(gx)]) - np.min(gx[~np.isnan(gx)])))
    gxn *= (gx - (
                (np.min(gx[~np.isnan(gx)]) + np.max(gx[~np.isnan(gx)])) / 2.))
    gxn[np.isnan(gx)] = np.nan
    return gxn


def strike_dip_vector(strike, dip):
    vec = np.zeros((len(strike), 3))
    s_r = np.deg2rad(strike)
    d_r = np.deg2rad((dip))
    vec[:, 0] = np.sin(d_r) * np.cos(s_r)
    vec[:, 1] = -np.sin(d_r) * np.sin(s_r)
    vec[:, 2] = np.cos(d_r)
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    return vec


def normal_vector_to_strike_and_dip(normal_vector):
    normal_vector /= np.linalg.norm(normal_vector, axis=1)[:, None]
    dip = np.rad2deg(np.arccos(normal_vector[:, 2]));
    strike = np.rad2deg(np.arctan2(normal_vector[:, 1], normal_vector[:,
                                                        0]))  # atan2(v2[1],v2[0])*rad2deg;
    return np.array([strike, dip]).T
