import numpy as np
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class BaseStructuredSupport:
    """ """

    def __init__(
        self,
        origin=np.zeros(3),
        nsteps=np.array([10, 10, 10]),
        step_vector=np.ones(3),
        rotation_xy=None,
    ):
        """

        Parameters
        ----------
        origin - 3d list or numpy array
        nsteps - 3d list or numpy array of ints
        step_vector - 3d list or numpy array of int
        """
        # the geometry in the mesh can be calculated from the
        # nsteps, step vector and origin
        # we use property decorators to update these when different parts of
        # the geometry need to change
        # inisialise the private attributes
        if np.any(step_vector == 0):
            logger.warning(f"Step vector {step_vector} has zero values")
        self._nsteps = np.array(nsteps, dtype=int) + 1
        self._step_vector = np.array(step_vector)
        self._origin = np.array(origin)
        self.supporttype = "Base"
        self._rotation_xy = np.zeros((3, 3))
        self._rotation_xy[0, 0] = 1
        self._rotation_xy[1, 1] = 1
        self._rotation_xy[2, 2] = 1
        self.rotation_xy = rotation_xy

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, nsteps):
        # if nsteps changes we need to change the step vector
        change_factor = nsteps / self.nsteps
        self._step_vector /= change_factor
        self._nsteps = nsteps

    @property
    def nsteps_cells(self):
        return self.nsteps - 1

    @property
    def rotation_xy(self):
        return self._rotation_xy

    @rotation_xy.setter
    def rotation_xy(self, rotation_xy):
        if rotation_xy is not None:
            self._rotation_xy[:, :] = 0
            self._rotation_xy[0, 0] = np.cos(np.deg2rad(rotation_xy))
            self._rotation_xy[0, 1] = -np.sin(np.deg2rad(rotation_xy))
            self._rotation_xy[1, 0] = np.sin(np.deg2rad(rotation_xy))
            self._rotation_xy[1, 1] = np.cos(np.deg2rad(rotation_xy))
            self._rotation_xy[2, 2] = 1.0  # make sure rotation around z vector

    @property
    def step_vector(self):
        return self._step_vector

    @step_vector.setter
    def step_vector(self, step_vector):
        change_factor = step_vector / self._step_vector
        newsteps = self._nsteps / change_factor
        self._nsteps = np.ceil(newsteps).astype(int)
        self._step_vector = step_vector

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        origin = np.array(origin)
        length = self.maximum - origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(int)
        self._origin = origin

    @property
    def maximum(self):
        return self.origin + self.nsteps_cells * self.step_vector

    @maximum.setter
    def maximum(self, maximum):
        """
        update the number of steps to fit new boundary
        """
        maximum = np.array(maximum, dtype=float)
        length = maximum - self.origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(int) + 1

    @property
    def n_nodes(self):
        return np.product(self.nsteps)

    @property
    def n_elements(self):
        return np.product(self.nsteps_cells)

    def __str__(self):
        return (
            "LoopStructural interpolation support:  {} \n"
            "Origin: {} {} {} \n"
            "Maximum: {} {} {} \n"
            "Step Vector: {} {} {} \n"
            "Number of Steps: {} {} {} \n"
            "Degrees of freedon {}".format(
                self.supporttype,
                self.origin[0],
                self.origin[1],
                self.origin[2],
                self.maximum[0],
                self.maximum[1],
                self.maximum[2],
                self.step_vector[0],
                self.step_vector[1],
                self.step_vector[2],
                self.nsteps[0],
                self.nsteps[1],
                self.nsteps[2],
                self.n_nodes,
            )
        )

    @property
    def nodes(self):
        max = self.origin + self.nsteps_cells * self.step_vector
        if np.any(np.isnan(self.nsteps)):
            raise ValueError("Cannot resize mesh nsteps is NaN")
        if np.any(np.isnan(self.origin)):
            raise ValueError("Cannot resize mesh origin is NaN")

        x = np.linspace(self.origin[0], max[0], self.nsteps[0])
        y = np.linspace(self.origin[1], max[1], self.nsteps[1])
        z = np.linspace(self.origin[2], max[2], self.nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        return np.array(
            [xx.flatten(order="F"), yy.flatten(order="F"), zz.flatten(order="F")]
        ).T

    def rotate(self, pos):
        """ """
        return np.einsum("ijk,ik->ij", self.rotation_xy[None, :, :], pos)

    def position_to_cell_index(self, pos):
        """Get the indexes (i,j,k) of a cell
        that a point is inside


        Parameters
        ----------
        pos : np.array
            Nx3 array of xyz locations

        Returns
        -------
        np.array, np.array, np.array
            i,j,k indexes of the cell that the point is in
        """

        pos = self.check_position(pos)

        ix = pos[:, 0] - self.origin[None, 0]
        iy = pos[:, 1] - self.origin[None, 1]
        iz = pos[:, 2] - self.origin[None, 2]
        ix = ix // self.step_vector[None, 0]
        iy = iy // self.step_vector[None, 1]
        iz = iz // self.step_vector[None, 2]
        return ix.astype(int), iy.astype(int), iz.astype(int)

    def position_to_cell_global_index(self, pos):
        ix, iy, iz = self.position_to_cell_index(pos)

    def inside(self, pos):
        # check whether point is inside box
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= pos[:, i] < self.maximum[None, i]
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

    def _global_indicies(self, indexes, nsteps):
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        if len(indexes.shape) == 1:
            raise ValueError("Cell indexes needs to be Nx3")
        if len(indexes.shape) == 2:
            if indexes.shape[1] != 3 and indexes.shape[0] == 3:
                indexes = indexes.swapaxes(0, 1)
            if indexes.shape[1] != 3:
                logger.error("Indexes shape {}".format(indexes.shape))
                raise ValueError("Cell indexes needs to be Nx3")
            return (
                indexes[:, 0]
                + nsteps[None, 0] * indexes[:, 1]
                + nsteps[None, 0] * nsteps[None, 1] * indexes[:, 2]
            )
        if len(indexes.shape) == 3:
            if indexes.shape[2] != 3 and indexes.shape[1] == 3:
                indexes = indexes.swapaxes(1, 2)
            if indexes.shape[2] != 3 and indexes.shape[0] == 3:
                indexes = indexes.swapaxes(0, 2)
            if indexes.shape[2] != 3:
                logger.error("Indexes shape {}".format(indexes.shape))
                raise ValueError("Cell indexes needs to be NxNx3")
            return (
                indexes[:, :, 0]
                + nsteps[None, None, 0] * indexes[:, :, 1]
                + nsteps[None, None, 0] * nsteps[None, None, 1] * indexes[:, :, 2]
            )

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

    def position_to_cell_corners(self, pos):

        inside = self.inside(pos)
        ix, iy, iz = self.position_to_cell_index(pos)
        cornersx, cornersy, cornersz = self.cell_corner_indexes(ix, iy, iz)
        globalidx = self.global_indicies(np.dstack([cornersx, cornersy, cornersz]).T)
        # if global index is not inside the support set to -1
        globalidx[~inside] = -1
        return globalidx, inside

    def position_to_cell_vertices(self, pos):
        """Get the vertices of the cell a point is in

        Parameters
        ----------
        pos : np.array
            Nx3 array of xyz locations

        Returns
        -------
        np.array((N,3),dtype=float), np.array(N,dtype=int)
            vertices, inside
        """
        gi, inside = self.position_to_cell_corners(pos)
        ci, cj, ck = self.global_index_to_node_index(gi)
        return self.node_indexes_to_position(ci, cj, ck), inside

    def node_indexes_to_position(self, xindex, yindex, zindex):

        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex
        z = self.origin[2] + self.step_vector[2] * zindex

        return x, y, z

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
        z_index = (
            global_index // self.nsteps_cells[0, None] // self.nsteps_cells[1, None]
        )
        return x_index, y_index, z_index

    def global_index_to_node_index(self, global_index):
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
        x_index = global_index % self.nsteps[0, None]
        y_index = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        z_index = global_index // self.nsteps[0, None] // self.nsteps[1, None]
        return x_index, y_index, z_index

    def global_node_indicies(self, indexes):
        """
        Convert from node indexes to global node index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._global_indicies(indexes, self.nsteps)

    def global_cell_indicies(self, indexes):
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._global_indicies(indexes, self.nsteps_cells)
