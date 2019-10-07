import numpy as np
from LoopStructural.supports.base_grid import BaseGrid

import logging
logger = logging.getLogger(__name__)

class StructuredGrid(BaseGrid):
    def __init__(self, nsteps=np.array([10, 10]), step_vector=np.ones(2), origin=np.zeros(2)):
        self.nsteps = nsteps
        self.step_vector = step_vector
        self.origin = origin
        self.n = self.nsteps[0] * self.nsteps[1]
        self.ndim = (2)  # make unmutable

    def position_to_cell_index(self, pos):
        """

        :param pos:
        :return:
        """
        pos = self.check_position(pos)
        # if type(pos) != np.array:
        #     return pos
        # inside = self.is_inside(pos)
        ix = pos[:, 0] - self.origin[None, 0]
        iy = pos[:, 1] - self.origin[None, 1]
        ix = ix.astype(int) // self.step_vector[None, 0]
        iy = iy.astype(int) // self.step_vector[None, 1]
        return ix.astype(int), iy.astype(int)

    def check_position(self, pos):
        """

        :param pos:
        :return:
        """
        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            print("Position array needs to be a list of points or a point")
            return False
        return pos

    def bilinear_interpolation(self, local_x, local_y):
        """

        :param local_x:
        :param local_y:
        :param local_z:
        :return:
        """

        return np.array([(1 - local_x) * (1 - local_y), local_x * (1 - local_y), (1 - local_x) * local_y,
                         local_x * local_y])

    def position_to_local_coordinates(self, pos):
        """

        :param pos:
        :return:
        """
        # TODO check if inside mesh
        # calculate local coordinates for positions
        local_x = ((pos[:, 0] - self.origin[None, 0]) % self.step_vector[None, 0]) / self.step_vector[None, 0]
        local_y = ((pos[:, 1] - self.origin[None, 1]) % self.step_vector[None, 1]) / self.step_vector[None, 1]

        return local_x, local_y

    def position_to_dof_coefs(self, pos):
        """

        :param pos:
        :return:
        """
        x_local, y_local = self.position_to_local_coordinates(pos)
        weights = self.bilinear_interpolation(x_local, y_local)
        return weights

    def global_indicies(self, indexes):
        """
        convert from xyz indexes to global index
        :param indexes:
        :return:
        """
        return indexes[:, 0] + self.nsteps[None, 0] * indexes[:, 1]

    def neighbour_global_indexes(self):
        """
        Find the neighbours= nodes for a node
        :param indexes: numpy array containing x,y,z indexes
        :return: n*9 array of global indexes
        """
        """return -1 for neighbours that would be on the overside of the border"""

        ii = []
        jj = []
        for i in range(1, self.nsteps[0] - 1):
            for j in range(1, self.nsteps[1] - 1):
                ii.append(i)
                jj.append(j)
        indexes = np.array([ii, jj])
        mask = np.array([[-1, 0, 1, -1, 0, 1, -1, 0, 1], [-1, -1, -1, 0, 0, 0, 1, 1, 1]])
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1, :, :]

    def cell_corner_indexes(self, x_cell_index, y_cell_index):
        """

        :param x_cell_index:
        :param y_cell_index:
        :param z_cell_index:
        :return:
        """
        xcorner = np.array([0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1])
        xcorners = x_cell_index[:, None] + xcorner[None, :]
        ycorners = y_cell_index[:, None] + ycorner[None, :]
        return xcorners, ycorners

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global index to grid index
        :param global_index:
        :return: x,y,z indices for the grid
        """
        x_index = global_index % self.nsteps[0, None]
        y_index = global_index // self.nsteps[0, None]
        return x_index, y_index

    def node_indexes_to_position(self, xindex, yindex):
        """
        convert array of xindexes and yindexes to positions
        :param xindex:
        :param yindex:
        :param zindex:

        :return: length Nx numpy array of coordinates
        """

        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex

        return x, y

    def position_to_cell_corners(self, pos):
        """
        find containing cell and get the corner coords
        :param pos:
        :return:
        """
        ix, iy = self.position_to_cell_index(pos)
        cornersx, cornersy = self.cell_corner_indexes(ix, iy)
        return cornersx + self.nsteps[0] * cornersy

    def calcul_T(self, pos):
        """
        Calculates the gradient matrix at location pos
        :param pos: numpy array of location Nx3
        :return: Nx3x4 matrix
        """
        xindex, yindex = self.position_to_cell_index(pos)
        cellx, celly = self.cell_corner_indexes(xindex, yindex)
        x, y = self.node_indexes_to_position(cellx, celly)
        T = np.zeros((pos.shape[0], 2, 4))

        div = self.step_vector[0] * self.step_vector[1]
        T[:, 0, 0] = -(y[:, 2] - pos[:, 1]) / div
        T[:, 0, 1] = (y[:, 3] - pos[:, 1]) / div
        T[:, 0, 2] = -(pos[:, 1] - y[:, 0]) / div
        T[:, 0, 3] = (pos[:, 1] - y[:, 1]) / div

        T[:, 1, 0] = -(x[:, 1] - pos[:, 0]) / div
        T[:, 1, 1] = -(pos[:, 0] - x[:, 0]) / div
        T[:, 1, 2] = (x[:, 3] - pos[:, 0]) / div
        T[:, 1, 3] = (pos[:, 0] - x[:, 2]) / div

        return T

    def assemble_inner(self, operator):
        """

        :param operator:
        :return:
        """
        A = []
        B = []
        cols = []
        global_indexes = self.neighbour_global_indexes()  # np.array([ii,jj]))
        a = np.tile(operator.flatten(), global_indexes.shape[1])
        r = np.array([np.arange(0, global_indexes.shape[1])])
        r += self.nc
        self.nc += global_indexes.shape[1]
        r = np.tile(r, (9, 1))
        r = r.T  # .flatten()
        r = r.flatten()
        A.extend(a.flatten().tolist())
        cols.extend(global_indexes.T.flatten().tolist())
        B.extend(np.zeros(global_indexes.shape[1]).tolist())
        return A,

    def assemble_borders(self, operator):
        """
        Add border masks to the interpolation matrix
        :return:
        """

        # TODO move this to the grid class
        # bottom and top borders
        mask = np.array([-1, 0, 1])
        jj = np.arange(1, self.grid.nsteps[0] - 1)
        ii = np.arange(1, self.grid.nsteps[1] - 1)
        top = (mask[None, :] + jj[:, None]).flatten() + self.grid.nsteps[0] * 0
        bottom = (mask[None, :] + jj[:, None]).flatten() + self.grid.nsteps[0] * (self.grid.nsteps[1] - 1)
        left = 0 + self.grid.nsteps[0] * (mask[None, :] + ii[:, None]).flatten()
        right = (self.grid.nsteps[0] - 1) + self.grid.nsteps[0] * (mask[None, :] + ii[:, None]).flatten()
        # add top and bottom to interpolation matrix
        a = np.tile(operator.flatten(), len(jj))
        r = np.array([np.arange(0, len(jj))])
        r += self.nc
        self.nc += len(jj)
        r = np.tile(r, (3, 1))
        r = r.T  # .flatten()
        r = r.flatten()
        self.A.extend(a.flatten().tolist())
        self.cols.extend(top.flatten().tolist())
        self.rows.extend(r.tolist())
        self.B.extend(np.zeros(len(jj)).tolist())
        r += len(jj)

        self.nc += len(jj)
        self.A.extend(a.flatten().tolist())
        self.cols.extend(bottom.flatten().tolist())
        self.rows.extend(r.tolist())
        self.B.extend(np.zeros(len(jj)).tolist())

        # add left and right to interpolation matrix
        a = np.tile(operator.flatten(), len(ii))
        r = np.array([np.arange(0, len(ii))])
        r += self.nc
        self.nc += len(ii)
        r = np.tile(r, (3, 1))
        r = r.T  # .flatten()
        r = r.flatten()
        self.A.extend(a.flatten().tolist())
        self.cols.extend(left.flatten().tolist())
        self.rows.extend(r.tolist())
        self.B.extend(np.zeros(len(ii)).tolist())

        r += len(ii)
        self.nc += len(ii)
        self.A.extend(a.flatten().tolist())
        self.cols.extend(right.flatten().tolist())
        self.rows.extend(r.tolist())
        self.B.extend(np.zeros(len(ii)).tolist())
#
#     self.origin = origin
#     self.shape = shape
#     self.step_vector=step_vector
#     self.n = shape[0]*shape[1]*shape[2]
#     self.n_cell_x = self.shape[0]-1
#     self.n_cell_y = self.shape[1] - 1
#     self.n_cell_z = self.shape[2] - 1
#     self.n_cell = self.n_cell_x*self.n_cell_y*self.n_cell_z
#     self.aa = True
#     if 'aa' in kwargs:
#         self.aa = kwargs['aa']
#     #TODO add in PCA for axis alignment
# def position_is_inside(self,pos):
#     #check if point is inside mesh
#
# def position_to_dof_coefs(self,x,y):
#
#     import numpy as np
#     class Grid2D:
#         """
#
#         """
#
#         def __init__(self, nsteps=np.array([10, 10]), step_vector=np.ones(2), origin=np.zeros(2)):
#             """
#
#             :param nsteps:
#             :param step_vector:
#             :param origin:
#             """
#             self.nsteps = nsteps
#             self.step_vector = step_vector
#             self.origin = origin
#             self.n = self.nsteps[0] * self.nsteps[1]
#             self.ndim = (2)  # make unmutable

