"""
Cartesian grid for fold interpolator

"""
import logging

import numpy as np

from .base_structured_3d_support import BaseStructuredSupport


from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class StructuredGrid(BaseStructuredSupport):
    """

    """
    def __init__(self,
                 origin=np.zeros(3),
                 nsteps=np.array([10, 10, 10]),
                 step_vector=np.ones(3),
                 name=None
                 ):
        """

        Parameters
        ----------
        origin - 3d list or numpy array
        nsteps - 3d list or numpy array of ints
        step_vector - 3d list or numpy array of int
        """
        BaseStructuredSupport.__init__(self,origin,nsteps,step_vector)
        self.regions = {}
        self.regions['everywhere'] = np.ones(self.n_nodes).astype(bool)
        self.name = name

    def barycentre(self):
        return self.cell_centres(np.arange(self.n_elements))


    def cell_centres(self, global_index):
        """get the centre of specified cells

        
        Parameters
        ----------
        global_index : array/list
            container of integer global indexes to cells

        Returns
        -------
        numpy array
            Nx3 array of cell centres
        """
        ix, iy, iz = self.global_index_to_cell_index(global_index)
        x = self.origin[None, 0] + self.step_vector[None, 0] * .5 + \
            self.step_vector[None, 0] * ix
        y = self.origin[None, 1] + self.step_vector[None, 1] * .5 + \
            self.step_vector[None, 1] * iy
        z = self.origin[None, 2] + self.step_vector[None, 2] * .5 + \
            self.step_vector[None, 2] * iz
        return np.array([x, y, z]).T





    def trilinear(self, x, y, z):
        """
        returns the trilinear interpolation for the local coordinates
        Parameters
        ----------
        x - double, array of doubles
        y - double, array of doubles
        z - double, array of doubles

        Returns
        -------
        array of interpolation coefficients

        """
        return np.array([(1 - x) * (1 - y) * (1 - z),
                         x * (1 - y) * (1 - z),
                         (1 - x) * y * (1 - z),
                         (1 - x) * (1 - y) * z,
                         x * (1 - y) * z,
                         (1 - x) * y * z,
                         x * y * (1 - z),
                         x * y * z])

    def position_to_local_coordinates(self, pos):
        """
        Convert from global to local coordinates within a cel
        Parameters
        ----------
        pos - array of positions inside

        Returns
        -------
        localx, localy, localz

        """
        # TODO check if inside mesh
        # pos = self.rotate(pos)
        # calculate local coordinates for positions
        local_x = ((pos[:, 0] - self.origin[None, 0]) % self.step_vector[
            None, 0]) / self.step_vector[None, 0]
        local_y = ((pos[:, 1] - self.origin[None, 1]) % self.step_vector[
            None, 1]) / self.step_vector[None, 1]
        local_z = ((pos[:, 2] - self.origin[None, 2]) % self.step_vector[
            None, 2]) / self.step_vector[None, 2]
        return local_x, local_y, local_z

    def position_to_dof_coefs(self, pos):
        """
        global posotion to interpolation coefficients
        Parameters
        ----------
        pos

        Returns
        -------

        """
        x_local, y_local, local_z = self.position_to_local_coordinates(pos)
        weights = self.trilinear(x_local, y_local, local_z)
        return weights

    def global_indicies(self, indexes):
        """
        xi, yi, zi to global index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        indexes = np.array(indexes).swapaxes(0, 2)
        return indexes[:, :, 0] + self.nsteps[None, None, 0] * indexes[:, :,
                                                               1] + \
               self.nsteps[None, None, 0] * self.nsteps[
                   None, None, 1] * indexes[:, :, 2]

    def neighbour_global_indexes(self, mask = None, **kwargs):
        """
        Get neighbour indexes

        Parameters
        ----------
        kwargs - indexes array specifying the cells to return neighbours

        Returns
        -------

        """
        indexes = None
        if "indexes" in kwargs:
            indexes = kwargs['indexes']
        if "indexes" not in kwargs:
            ii = []
            jj = []
            kk = []
            for i in range(1, self.nsteps[0] - 1):
                for j in range(1, self.nsteps[1] - 1):
                    for k in range(1, self.nsteps[2] - 1):
                        kk.append(k)
                        ii.append(i)
                        jj.append(j)
            indexes = np.array([ii, jj, kk])
        # indexes = np.array(indexes).T
        if indexes.ndim != 2:
            print(indexes.ndim)
            return
        # determine which neighbours to return default is diagonals included.
        if mask is None:
            mask = np.array([
                [-1, 0, 1, -1, 0, 1, -1, 0, 1,
                 -1, 0, 1, -1, 0, 1, -1, 0, 1,
                 -1, 0, 1, -1, 0, 1, -1, 0, 1],
                [-1, -1, -1, 0, 0, 0, 1, 1, 1,
                 -1, -1, -1, 0, 0, 0, 1, 1, 1,
                 -1, -1, -1, 0, 0, 0, 1, 1, 1],
                [-1, -1, -1, -1, -1, -1, -1, -1, -1,
                 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 1, 1, 1, 1, 1]
            ])
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return(neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1,
                                                                  :, :] + \
               self.nsteps[0, None, None] * self.nsteps[
                   1, None, None] * neighbours[2, :, :]).astype(np.int64)

    def cell_corner_indexes(self, x_cell_index, y_cell_index, z_cell_index):
        """
        Returns the indexes of the corners of a cell given its location xi,
        yi, zi

        Parameters
        ----------
        x_cell_index
        y_cell_index
        z_cell_index

        Returns
        -------

        """
        xcorner = np.array([0, 1, 0, 0, 1, 0, 1, 1])
        ycorner = np.array([0, 0, 1, 0, 0, 1, 1, 1])
        zcorner = np.array([0, 0, 0, 1, 1, 1, 0, 1])
        xcorners = x_cell_index[:, None] + xcorner[None, :]
        ycorners = y_cell_index[:, None] + ycorner[None, :]
        zcorners = z_cell_index[:, None] + zcorner[None, :]
        return xcorners, ycorners, zcorners

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        # determine the ijk indices for the global index.
        # remainder when dividing by nx = i
        # remained when dividing modulus of nx by ny is j

        x_index = global_index % self.nsteps_cells[0, None]
        y_index = global_index // self.nsteps_cells[0, None] % \
                  self.nsteps_cells[1, None]
        z_index = global_index // self.nsteps_cells[0, None] // \
                  self.nsteps_cells[1, None]
        return x_index, y_index, z_index

    def node_indexes_to_position(self, xindex, yindex, zindex):

        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex
        z = self.origin[2] + self.step_vector[2] * zindex

        return x, y, z

    def position_to_cell_corners(self, pos):

        inside = self.inside(pos)
        ix, iy, iz = self.position_to_cell_index(pos)
        cornersx, cornersy, cornersz = self.cell_corner_indexes(ix, iy, iz)
        globalidx = self.global_indicies(
            np.dstack([cornersx, cornersy, cornersz]).T)
        # if global index is not inside the support set to -1
        globalidx[~inside] = -1
        return globalidx, inside

    def evaluate_value(self, evaluation_points, property_array):
        """
        Evaluate the value of of the property at the locations.
        Trilinear interpolation dot corner values

        Parameters
        ----------
        evaluation_points np array of locations
        property_name string of property name

        Returns
        -------

        """
        if property_array.shape[0] != self.n_nodes:
            logger.error("Property array does not match grid")
            raise BaseException
        idc, inside = self.position_to_cell_corners(evaluation_points)
        v = np.zeros(idc.shape)
        v[:, :] = np.nan

        v[inside, :] = self.position_to_dof_coefs(
            evaluation_points[inside, :]).T
        
        v[inside, :] *= property_array[idc[inside, :]]
        
        return np.sum(v, axis=1)

    def evaluate_gradient(self, evaluation_points, property_array):
        if property_array.shape[0] != self.n_nodes:
            logger.error("Property array does not match grid")
            raise BaseException
        idc, inside = self.position_to_cell_corners(evaluation_points)
        T = np.zeros((idc.shape[0], 3, 8))
        T[inside, :, :] = self.calcul_T(evaluation_points[inside, :])
        # indices = np.array([self.position_to_cell_index(evaluation_points)])
        # idc = self.global_indicies(indices.swapaxes(0,1))
        # print(idc)
        T[inside, 0, :] *= property_array[idc[inside, :]]
        T[inside, 1, :] *= property_array[idc[inside, :]]
        T[inside, 2, :] *= property_array[idc[inside, :]]
        return np.array(
            [np.sum(T[:, 0, :], axis=1), np.sum(T[:, 1, :], axis=1) ,
             np.sum(T[:, 2, :], axis=1) ]).T

    def calcul_T(self, pos):
        """
        Calculates the gradient matrix at location pos
        :param pos: numpy array of location Nx3
        :return: Nx3x4 matrix
        """
        #   6_ _ _ _ 8
        #   /|    /|
        # 4 /_|  5/ |
        # | 2|_ _|_| 7
        # | /    | /
        # |/_ _ _|/
        # 0      1
        #
        # xindex, yindex, zindex = self.position_to_cell_index(pos)
        # cellx, celly, cellz = self.cell_corner_indexes(xindex, yindex,zindex)
        # x, y, z = self.node_indexes_to_position(cellx, celly, cellz)
        T = np.zeros((pos.shape[0], 3, 8))
        x, y, z = self.position_to_local_coordinates(pos)
    
        T[:, 0, 0] = (1 - z) * (y- 1)  # v000
        T[:, 0, 1] = (1 - y) * (1 - z)  # (y[:, 3] - pos[:, 1]) / div
        T[:, 0, 2] = -y * (1 - z)  # (pos[:, 1] - y[:, 0]) / div
        T[:, 0, 3] = -(1 - y) * z  # (pos[:, 1] - y[:, 1]) / div
        T[:, 0, 4] = (1 - y) * z
        T[:, 0, 5] = - y * z
        T[:, 0, 6] = y * (1 - z)
        T[:, 0, 7] = y * z

        T[:, 1, 0] =  (x - 1) * (1 - z)
        T[:, 1, 1] = - x * (1 - z)
        T[:, 1, 2] = (1 - x) * (1 - z)
        T[:, 1, 3] = -(1 - x) * z
        T[:, 1, 4] = -x * z
        T[:, 1, 5] = (1 - x) * z
        T[:, 1, 6] = x * (1 - z)
        T[:, 1, 7] = x * z

        T[:, 2, 0] = -(1 - x) * (1 - y)
        T[:, 2, 1] = - x * (1 - y)
        T[:, 2, 2] = - (1 - x) * y
        T[:, 2, 3] = (1 - x) * (1 - y)
        T[:, 2, 4] = x * (1 - y)
        T[:, 2, 5] = (1 - x) * y
        T[:, 2, 6] = - x * y
        T[:, 2, 7] = x * y
        T/=self.step_vector[0]

        return T 

