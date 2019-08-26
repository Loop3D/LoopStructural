import numpy as np
from .base_grid import BaseGrid
from skimage.measure import marching_cubes_lewiner as marching_cubes
class StructuredGrid:
    def __init__(self,
                 nsteps=np.array([10, 10, 10]),
                 step_vector=np.ones(3),
                 origin=np.zeros(3)):
        """

        :param nsteps:
        :param step_vector:
        :param origin:
        """
        self.nsteps = np.array(nsteps)
        self.step_vector = np.array(step_vector)
        self.origin = np.array(origin)
        self.n = self.nsteps[0] * self.nsteps[1] * self.nsteps[2]
        self.dim = 3
        self.n_cell_x = self.nsteps[0] - 1
        self.n_cell_y = self.nsteps[1] - 1
        self.n_cell_z = self.nsteps[2] - 1
        self.properties = {}
        self.n_cell = self.n_cell_x * self.n_cell_y * self.n_cell_z
        # BaseGrid.__init__(np.zeros((4,2)))


    def update_property(self, propertyname, values):
        if values.shape[0] == self.n:
            self.properties[propertyname] = values
        if values.shape[0] == self.n_cell:
            self.cell_properties[propertyname] = values

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
        iz = pos[:, 2] - self.origin[None, 2]
        ix = ix.astype(int) // self.step_vector[None, 0]
        iy = iy.astype(int) // self.step_vector[None, 1]
        iz = iz.astype(int) // self.step_vector[None, 2]
        return ix.astype(int), iy.astype(int), iz.astype(int)

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

    def trilinear(self, x, y, z):
        """
        trilinear interpolation for a unit cube
        :param x:
        :param y:
        :param z:
        :return:
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

        :param pos:
        :return:
        """
        # TODO check if inside mesh
        # calculate local coordinates for positions
        local_x = ((pos[:, 0] - self.origin[None, 0]) % self.step_vector[None, 0]) / self.step_vector[None, 0]
        local_y = ((pos[:, 1] - self.origin[None, 1]) % self.step_vector[None, 1]) / self.step_vector[None, 1]
        local_z = ((pos[:, 2] - self.origin[None, 2]) % self.step_vector[None, 2]) / self.step_vector[None, 2]

        return local_x, local_y, local_z

    def position_to_dof_coefs(self, pos):
        """

        :param pos:
        :return:
        """
        x_local, y_local, local_z = self.position_to_local_coordinates(pos)
        weights = self.trilinear(x_local, y_local, local_z)
        return weights

    def global_indicies(self, indexes):
        """
        convert from xyz indexes to global index
        :param indexes:
        :return:
        """
        indexes = np.array(indexes).swapaxes(0,2)

        return indexes[:,:, 0] + self.nsteps[None,None, 0] * indexes[:,:, 1] + \
               self.nsteps[None,None, 0] * self.nsteps[None,None, 1] * indexes[:, :, 2]

    def neighbour_global_indexes(self, **kwargs):
        """
        Find the neighbours= nodes for a node
        :param indexes: numpy array containing x,y,z indexes
        :return: n*9 array of global indexes
        """
        """return -1 for neighbours that would be on the overside of the border"""
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
        # neighbours = neighbours.swapaxes(0,1)
        # #calculate a mask to tell which neighbours are out of bounds
        # mask = neighbours[:,0,:]>-1
        # mask = np.logical_and(mask, neighbours[:,1,:]>-1)
        # mask = np.logical_and(mask, neighbours[:,2,:]>-1)
        # mask = np.logical_and(mask, neighbours[:,0,:]<self.nsteps[0])
        # mask = np.logical_and(mask, neighbours[:,1,:]<self.nsteps[1])
        # mask = np.logical_and(mask, neighbours[:,2,:]<self.nsteps[2])
        return neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1, :, :] + \
               self.nsteps[0, None, None] * self.nsteps[1, None, None] * neighbours[2, :, :]

    def cell_corner_indexes(self, x_cell_index, y_cell_index, z_cell_index):
        """

        :param x_cell_index:
        :param y_cell_index:
        :param z_cell_index:
        :return:
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
        Convert from global index to grid index
        :param global_index:
        :return: x,y,z indices for the grid
        """
        x_index = global_index % self.nsteps[0, None]
        y_index = global_index // self.nsteps[1, None] % self.nsteps[0, None]
        z_index = global_index // self.nsteps[0, None] // self.nsteps[1, None]
        return x_index, y_index, z_index

    def node_indexes_to_position(self, xindex, yindex, zindex):
        """
        convert array of xindexes and yindexes to positions
        :param xindex:
        :param yindex:
        :param zindex:

        :return: length Nx numpy array of coordinates
        """

        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex
        z = self.origin[2] + self.step_vector[2] * zindex

        return x, y, z

    def position_to_cell_corners(self, pos):
        """
        find containing cell and get the corner coords
        :param pos:
        :return:
        """
        ix, iy, iz = self.position_to_cell_index(pos)
        cornersx, cornersy, cornersz = self.cell_corner_indexes(ix, iy, iz)
        globalidx = self.global_indicies(np.dstack([cornersx, cornersy, cornersz]).T)
        return globalidx

    def evaluate_value(self, evaluation_points, property_name):
        corners = self.position_to_cell_corners(evaluation_points)
        dof = self.position_to_dof_coefs(evaluation_points)

    def calcul_T(self, pos):
        """
        Calculates the gradient matrix at location pos
        :param pos: numpy array of location Nx3
        :return: Nx3x4 matrix
        """
        #   6_ _ _ _ 8
        #   /|    /|
        #4 /_|  5/ |
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
        # div = self.step_vector[0] * self.step_vector[1] * self.step_vector[2]

        T[:, 0, 0] = -(1 - y) * (1 - z)  #v000
        T[:, 0, 1] = (1 - y) * (1 - z)  # (y[:, 3] - pos[:, 1]) / div
        T[:, 0, 2] = -y * (1 - z)  # (pos[:, 1] - y[:, 0]) / div
        T[:, 0, 3] = -(1 - y) * z  # (pos[:, 1] - y[:, 1]) / div
        T[:, 0, 4] = (1 - y) * z
        T[:, 0, 5] = - y * z
        T[:, 0, 6] = y * (1 - z)
        T[:, 0, 7] = y * z

        T[:, 1, 0] = - (1 - x) * (1 - z)
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

        xcorner = np.array([0, 1, 1, 0, 0, 1, 1, 0])
        ycorner = np.array([0, 0, 1, 1, 0, 0, 1, 1])
        zcorner = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        return T

    def slice(self,propertyname, isovalue):
        verts, faces, normals, values = marching_cubes(
            self.properties[propertyname].reshape(self.nsteps,order='F'),
            isovalue,
            spacing=self.step_vector)
        return faces, verts