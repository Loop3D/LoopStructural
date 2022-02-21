"""
Cartesian grid for fold interpolator

"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class StructuredGrid2D:
    """ """

    def __init__(
        self,
        origin=np.zeros(2),
        nsteps=np.array([10, 10]),
        step_vector=np.ones(2),
    ):
        """

        Parameters
        ----------
        origin - 2d list or numpy array
        nsteps - 2d list or numpy array of ints
        step_vector - 2d list or numpy array of int
        """

        self.nsteps = np.array(nsteps)
        self.step_vector = np.array(step_vector)
        self.origin = np.array(origin)
        self.maximum = origin + self.nsteps * self.step_vector

        self.n_nodes = self.nsteps[0] * self.nsteps[1]
        self.dim = 2
        self.nsteps_cells = self.nsteps - 1
        self.n_cell_x = self.nsteps[0] - 1
        self.n_cell_y = self.nsteps[1] - 1
        self.properties = {}
        self.n_elements = self.n_cell_x * self.n_cell_y

        # calculate the node positions using numpy (this should probably not
        # be stored as it defeats
        # the purpose of a structured grid

        # self.barycentre = self.cell_centres(np.arange(self.n_elements))

        self.regions = {}
        self.regions["everywhere"] = np.ones(self.n_nodes).astype(bool)

    @property
    def nodes(self):
        max = self.origin + self.nsteps_cells * self.step_vector
        x = np.linspace(self.origin[0], max[0], self.nsteps[0])
        y = np.linspace(self.origin[1], max[1], self.nsteps[1])
        xx, yy = np.meshgrid(x, y, indexing="ij")
        return np.array([xx.flatten(order="F"), yy.flatten(order="F")]).T

    @property
    def barycentre(self):
        return self.cell_centres(np.arange(self.n_elements))

    # @property
    # def barycentre(self):
    #     return self.cell_centres(np.arange(self.n_elements))

    def print_geometry(self):
        print("Origin: %f %f %f" % (self.origin[0], self.origin[1], self.origin[2]))
        print(
            "Cell size: %f %f %f"
            % (self.step_vector[0], self.step_vector[1], self.step_vector[2])
        )
        max = self.origin + self.nsteps_cells * self.step_vector
        print("Max extent: %f %f %f" % (max[0], max[1], max[2]))

    def update_property(self, propertyname, values):
        """[summary]

        [extended_summary]

        Parameters
        ----------
        propertyname : [type]
            [description]
        values : [type]
            [description]
        """
        if values.shape[0] == self.n_nodes:
            self.properties[propertyname] = values
        if values.shape[0] == self.n_elements:
            self.cell_properties[propertyname] = values

    def cell_centres(self, global_index):
        """[summary]

        [extended_summary]

        Parameters
        ----------
        global_index : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        ix, iy = self.global_index_to_cell_index(global_index)
        x = (
            self.origin[None, 0]
            + self.step_vector[None, 0] * 0.5
            + self.step_vector[None, 0] * ix
        )
        y = (
            self.origin[None, 1]
            + self.step_vector[None, 1] * 0.5
            + self.step_vector[None, 1] * iy
        )
        return np.array([x, y]).T

    def position_to_cell_index(self, pos):
        """[summary]

        [extended_summary]

        Parameters
        ----------
        pos : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """

        pos = self.check_position(pos)

        ix = pos[:, 0] - self.origin[None, 0]
        iy = pos[:, 1] - self.origin[None, 1]
        ix = ix // self.step_vector[None, 0]
        iy = iy // self.step_vector[None, 1]
        return ix.astype(int), iy.astype(int)

    def inside(self, pos):

        # check whether point is inside box
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(self.dim):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i]
                + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside

    def check_position(self, pos):
        """[summary]

        [extended_summary]

        Parameters
        ----------
        pos : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """

        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            print("Position array needs to be a list of points or a point")
            return False
        return pos

    def bilinear(self, x, y):
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
        return np.array([(1 - x) * (1 - y), x * (1 - y), (1 - x) * y, x * y])

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

        # calculate local coordinates for positions
        local_x = (
            (pos[:, 0] - self.origin[None, 0]) % self.step_vector[None, 0]
        ) / self.step_vector[None, 0]
        local_y = (
            (pos[:, 1] - self.origin[None, 1]) % self.step_vector[None, 1]
        ) / self.step_vector[None, 1]

        return local_x, local_y

    def position_to_dof_coefs(self, pos):
        """
        global posotion to interpolation coefficients
        Parameters
        ----------
        pos

        Returns
        -------

        """
        x_local, y_local = self.position_to_local_coordinates(pos)
        weights = self.bilinear(x_local, y_local)
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
        return indexes[:, :, 0] + self.nsteps[None, None, 0] * indexes[:, :, 1]

    def neighbour_global_indexes(self, mask=None, **kwargs):
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
            indexes = kwargs["indexes"]
        if "indexes" not in kwargs:
            ii = []
            jj = []
            for i in range(1, self.nsteps[0] - 1):
                for j in range(1, self.nsteps[1] - 1):
                    ii.append(i)
                    jj.append(j)
            indexes = np.array([ii, jj])
        # indexes = np.array(indexes).T
        if indexes.ndim != 2:
            print(indexes.ndim)
            return
        # determine which neighbours to return default is diagonals included.
        if mask is None:
            mask = np.array(
                [[-1, 0, 1, -1, 0, 1, -1, 0, 1], [1, 1, 1, 0, 0, 0, -1, -1, -1]]
            )
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return (
            neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1, :, :]
        ).astype(np.int64)

    def cell_corner_indexes(self, x_cell_index, y_cell_index):
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
        xcorner = np.array([0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1])
        xcorners = x_cell_index[:, None] + xcorner[None, :]
        ycorners = y_cell_index[:, None] + ycorner[None, :]
        return xcorners, ycorners

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
        y_index = (
            global_index // self.nsteps_cells[0, None] % self.nsteps_cells[1, None]
        )
        return x_index, y_index

    def node_indexes_to_position(self, xindex, yindex):

        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex

        return x, y

    def position_to_cell_corners(self, pos):

        inside = self.inside(pos)
        ix, iy = self.position_to_cell_index(pos)
        cornersx, cornersy = self.cell_corner_indexes(ix, iy)
        globalidx = self.global_indicies(np.dstack([cornersx, cornersy]).T)
        # if global index is not inside the support set to -1
        globalidx[~inside] = -1
        return globalidx, inside

    def evaluate_value(self, evaluation_points, property_name):
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
        idc, inside = self.position_to_cell_corners(evaluation_points)
        v = np.zeros(idc.shape)
        v[:, :] = np.nan

        v[inside, :] = self.position_to_dof_coefs(evaluation_points[inside, :]).T
        v[inside, :] *= self.properties[property_name][idc[inside, :]]
        return np.sum(v, axis=1)

    def evaluate_gradient(self, evaluation_points, property_name):
        idc, inside = self.position_to_cell_corners(evaluation_points)
        T = np.zeros((idc.shape[0], 2, 4))
        T[inside, :, :] = self.calcul_T(evaluation_points[inside, :])
        # indices = np.array([self.position_to_cell_index(evaluation_points)])
        # idc = self.global_indicies(indices.swapaxes(0,1))
        # print(idc)
        T[inside, 0, :] *= self.properties[property_name][idc[inside, :]]
        T[inside, 1, :] *= self.properties[property_name][idc[inside, :]]
        # T[inside, 2, :] *= self.properties[property_name][idc[inside, :]]
        return np.array([np.sum(T[:, 0, :], axis=1), np.sum(T[:, 1, :], axis=1)]).T

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
        T = np.zeros((pos.shape[0], 2, 4))
        x, y = self.position_to_local_coordinates(pos)

        xindex, yindex = self.position_to_cell_index(pos)
        cellx, celly = self.cell_corner_indexes(xindex, yindex)
        x, y = self.node_indexes_to_position(cellx, celly)
        # div = self.step_vector[0] * self.step_vector[1] * self.step_vector[2]
        div = self.step_vector[0] * self.step_vector[1]

        T[:, 0, 0] = -(y[:, 2] - pos[:, 1]) / div
        T[:, 0, 1] = (y[:, 3] - pos[:, 1]) / div
        T[:, 0, 2] = -(pos[:, 1] - y[:, 0]) / div
        T[:, 0, 3] = (pos[:, 1] - y[:, 1]) / div

        T[:, 1, 0] = -(x[:, 1] - pos[:, 0]) / div
        T[:, 1, 1] = -(pos[:, 0] - x[:, 0]) / div
        T[:, 1, 2] = (x[:, 3] - pos[:, 0]) / div
        T[:, 1, 3] = (pos[:, 0] - x[:, 2]) / div

        # T[:, 0, 0] = -(1 - y)  # v000
        # T[:, 0, 1] = (1 - y)   # (y[:, 3] - pos[:, 1]) / div
        # T[:, 0, 2] = -y   # (pos[:, 1] - y[:, 0]) / div
        # T[:, 0, 3] = -(1 - y)  # (pos[:, 1] - y[:, 1]) / div

        # T[:, 1, 0] = - (1 - x)
        # T[:, 1, 1] = - x
        # T[:, 1, 2] = (1 - x)
        # T[:, 1, 3] = -(1 - x)

        return T
