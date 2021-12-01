import logging

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


def get_data_axis_aligned_bounding_box(xyz, buffer):
    """

    Parameters
    ----------
    xyz
    buffer

    Returns
    -------

    """
    minx = np.min(xyz[:, 0])
    maxx = np.max(xyz[:, 0])
    miny = np.min(xyz[:, 1])
    maxy = np.max(xyz[:, 1])
    minz = np.min(xyz[:, 2])
    maxz = np.max(xyz[:, 2])

    xlen = maxx - minx
    ylen = maxy - miny
    zlen = maxz - minz
    length = np.max([xlen, ylen, zlen])
    minx -= length * buffer
    maxx += length * buffer

    miny -= length * buffer
    maxy += length * buffer

    minz -= length * buffer
    maxz += length * buffer

    bb = np.array([[minx, miny, minz],
                   [maxx, maxy, maxz]
                   ])

    def region(xyz):
        """

        Parameters
        ----------
        xyz

        Returns
        -------

        """
        # print(xyz)
        # print(bb)
        b = np.ones(xyz.shape[0]).astype(bool)
        b = np.logical_and(b, xyz[:, 0] > minx)
        b = np.logical_and(b, xyz[:, 0] < maxx)
        b = np.logical_and(b, xyz[:, 1] > miny)
        b = np.logical_and(b, xyz[:, 1] < maxy)
        # b = np.logical_and(b,xyz[:,2]>minz)
        # b = np.logical_and(b,xyz[:,2]<maxz)
        return b

    return bb, region

def get_data_bounding_box_map(xyz, buffer):
    """

    Parameters
    ----------
    xyz
    buffer

    Returns
    -------

    """
    # find the aligned coordinates box using pca
    modelpca = PCA(n_components=3)
    modelpca.fit(xyz)
    # transform the data to this new coordinate then find extents
    transformed_xyz = modelpca.transform(xyz)
    minx = np.min(xyz[:, 0])
    maxx = np.max(xyz[:, 0])
    miny = np.min(xyz[:, 1])
    maxy = np.max(xyz[:, 1])
    minz = np.min(xyz[:, 2])
    maxz = np.max(xyz[:, 2])

    xlen = maxx - minx
    ylen = maxy - miny
    zlen = maxz - minz
    length = np.max([xlen, ylen, zlen])
    minx -= buffer
    maxx += buffer

    miny -= buffer
    maxy += buffer

    minz -= buffer
    maxz += buffer

    bb = np.array([[minx, miny, minz],
                   [maxx, maxy, maxz]
                   ])

    def region(xyz):
        b = np.ones(xyz.shape[0]).astype(bool)
        b = np.logical_and(b, xyz[:, 0] > minx)
        b = np.logical_and(b, xyz[:, 0] < maxx)
        b = np.logical_and(b, xyz[:, 1] > miny)
        b = np.logical_and(b, xyz[:, 1] < maxy)

        return b

    return bb, region
def get_data_bounding_box(xyz, buffer):
    """

    Parameters
    ----------
    xyz
    buffer

    Returns
    -------

    """
    # find the aligned coordinates box using pca
    modelpca = PCA(n_components=3)
    modelpca.fit(xyz)
    # transform the data to this new coordinate then find extents
    transformed_xyz = modelpca.transform(xyz)
    minx = np.min(xyz[:, 0])
    maxx = np.max(xyz[:, 0])
    miny = np.min(xyz[:, 1])
    maxy = np.max(xyz[:, 1])
    minz = np.min(xyz[:, 2])
    maxz = np.max(xyz[:, 2])

    xlen = maxx - minx
    ylen = maxy - miny
    zlen = maxz - minz
    length = np.max([xlen,ylen,zlen])
    minx -= length * buffer
    maxx += length * buffer

    miny -= length * buffer
    maxy += length * buffer

    minz -= length * buffer
    maxz += length * buffer

    bb = np.array([[minx, miny, minz],
                   [maxx, maxy, maxz]
                   ])

    def region(xyz):
        b = np.ones(xyz.shape[0]).astype(bool)
        b = np.logical_and(b, xyz[:, 0] > minx)
        b = np.logical_and(b, xyz[:, 0] < maxx)
        b = np.logical_and(b, xyz[:, 1] > miny)
        b = np.logical_and(b, xyz[:, 1] < maxy)
        b = np.logical_and(b, xyz[:, 2] > minz)
        b = np.logical_and(b, xyz[:, 2] < maxz)
        return b

    return bb, region


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

def create_box(bounding_box, nsteps):
    tri, xx, yy = create_surface(bounding_box[0:2, :], nsteps[0:2])

    zz = np.zeros(xx.shape)
    zz[:] = bounding_box[1, 2]

    tri = np.vstack([tri, tri + np.max(tri) + 1])
    xx = np.hstack([xx, xx])
    yy = np.hstack([yy, yy])

    z = np.zeros(zz.shape)
    z[:] = bounding_box[0, 2]
    zz = np.hstack([zz, z])
    # y faces
    t, x, z = create_surface(bounding_box[:, [0, 2]], nsteps[[0, 2]])
    tri = np.vstack([tri, t + np.max(tri) + 1])

    y = np.zeros(x.shape)
    y[:] = bounding_box[0, 1]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])

    tri = np.vstack([tri, t + np.max(tri) + 1])
    y[:] = bounding_box[1, 1]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])

    # x faces
    t, y, z = create_surface(bounding_box[:, [1, 2]], nsteps[[1, 2]])
    tri = np.vstack([tri, t + np.max(tri) + 1])
    x = np.zeros(y.shape)
    x[:] = bounding_box[0, 0]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])

    tri = np.vstack([tri, t + np.max(tri) + 1])
    x[:] = bounding_box[1, 0]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])

    points = np.zeros((len(xx), 3))  #
    points[:, 0] = xx
    points[:, 1] = yy
    points[:, 2] = zz
    return points, tri

def get_vectors(normal):
    length =  np.linalg.norm(normal,axis=1)[:,None]
    normal /= length#np.linalg.norm(normal,axis=1)[:,None]
    strikedip = normal_vector_to_strike_and_dip(normal)
    strike_vec = get_strike_vector(strikedip[:, 0])
    strike_vec /= np.linalg.norm(strike_vec,axis=0)[None,:]
    dip_vec = np.cross(strike_vec, normal, axisa=0, axisb=1).T  # (strikedip[:, 0], strikedip[:, 1])
    dip_vec /= np.linalg.norm(dip_vec,axis=0)[None,:]
    return strike_vec*length.T, dip_vec*length.T


def get_strike_vector(strike):
    """

    Parameters
    ----------
    strike

    Returns
    -------

    """

    v = np.array([np.sin(np.deg2rad(-strike)),
                  -np.cos(np.deg2rad(-strike)),
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
    dip = np.rad2deg(np.arccos(normal_vector[:, 2]))
    strike = -np.rad2deg(np.arctan2(normal_vector[:, 1], normal_vector[:,
                                                         0]))  # atan2(v2[1],v2[0])*rad2deg;

    return np.array([strike, dip]).T


def xyz_names():
    return ['X', 'Y', 'Z']


def normal_vec_names():
    return ['nx', 'ny', 'nz']


def tangent_vec_names():
    return ['tx', 'ty', 'tz']


def gradient_vec_names():
    return ['gx', 'gy', 'gz']


def weight_name():
    return ['w']


def val_name():
    return ['val']

def coord_name():
    return ['coord']

def interface_name():
    return ['interface']


def feature_name():
    return ['feature_name']

def all_heading():
    return xyz_names() + normal_vec_names() + tangent_vec_names() + \
           gradient_vec_names() + weight_name() + val_name() + coord_name() + feature_name() + interface_name()
def empty_dataframe():
    return pd.DataFrame(columns=[all_heading()])