from LoopStructural.utils.exceptions import LoopException
from abc import abstractmethod
import numpy as np
from LoopStructural.utils import getLogger
from . import SupportType
from typing import Tuple

logger = getLogger(__name__)

from ._base_support import BaseSupport


class BaseStructuredSupport(BaseSupport):
    """ """

    dimension = 3

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
        # cast to numpy array, to allow list like input
        origin = np.array(origin)
        nsteps = np.array(nsteps)
        step_vector = np.array(step_vector)

        self.type = SupportType.BaseStructured
        if np.any(step_vector == 0):
            logger.warning(f"Step vector {step_vector} has zero values")
        if np.any(nsteps == 0):
            raise LoopException("nsteps cannot be zero")
        if np.any(nsteps < 0):
            raise LoopException("nsteps cannot be negative")
        # if np.any(nsteps < 3):
        #     raise LoopException(
        #         "step vector cannot be less than 3. Try increasing the resolution of the interpolator"
        #     )
        self._nsteps = np.array(nsteps, dtype=int) + 1
        self._step_vector = np.array(step_vector)
        self._origin = np.array(origin)
        self.supporttype = "Base"
        self._rotation_xy = np.zeros((3, 3))
        self._rotation_xy[0, 0] = 1
        self._rotation_xy[1, 1] = 1
        self._rotation_xy[2, 2] = 1
        self.rotation_xy = rotation_xy
        self.interpolator = None

    @property
    def volume(self):
        return np.prod(self.maximum - self.origin)

    def set_nelements(self, nelements) -> int:
        box_vol = self.volume
        ele_vol = box_vol / nelements
        # calculate the step vector of a regular cube
        step_vector = np.zeros(3)

        step_vector[:] = ele_vol ** (1.0 / 3.0)

        # number of steps is the length of the box / step vector
        nsteps = np.ceil((self.maximum - self.origin) / step_vector).astype(int)
        self.nsteps = nsteps
        return self.n_elements

    def to_dict(self):
        return {
            "origin": self.origin,
            "nsteps": self.nsteps,
            "step_vector": self.step_vector,
            "rotation_xy": self.rotation_xy,
        }

    @abstractmethod
    def onGeometryChange(self):
        """Function to be called when the geometry of the support changes"""
        pass

    def associateInterpolator(self, interpolator):
        self.interpolator = interpolator

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, nsteps):
        # if nsteps changes we need to change the step vector
        change_factor = nsteps / self.nsteps
        self._step_vector /= change_factor
        self._nsteps = nsteps
        self.onGeometryChange()

    @property
    def nsteps_cells(self):
        return self.nsteps - 1

    @property
    def rotation_xy(self):
        return self._rotation_xy

    @rotation_xy.setter
    def rotation_xy(self, rotation_xy):
        if rotation_xy is None:
            return
        if isinstance(rotation_xy, (float, int)):
            rotation_xy = np.array(
                [
                    [
                        np.cos(np.deg2rad(rotation_xy)),
                        -np.sin(np.deg2rad(rotation_xy)),
                        0,
                    ],
                    [
                        np.sin(np.deg2rad(rotation_xy)),
                        np.cos(np.deg2rad(rotation_xy)),
                        0,
                    ],
                    [0, 0, 1],
                ]
            )
        rotation_xy = np.array(rotation_xy)
        if rotation_xy.shape != (3, 3):
            raise ValueError("Rotation matrix should be 3x3, not {}".format(rotation_xy.shape))
        self._rotation_xy = rotation_xy

    @property
    def step_vector(self):
        return self._step_vector

    @step_vector.setter
    def step_vector(self, step_vector):
        change_factor = step_vector / self._step_vector
        newsteps = self._nsteps / change_factor
        self._nsteps = np.ceil(newsteps).astype(int)
        self._step_vector = step_vector
        self.onGeometryChange()

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        origin = np.array(origin)
        length = self.maximum - origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(np.int64)
        self._nsteps[self._nsteps == 0] = (
            3  # need to have a minimum of 3 elements to apply the finite difference mask
        )
        if np.any(~(self._nsteps > 0)):
            logger.error(
                f"Cannot resize the interpolation support. The proposed number of steps is {self._nsteps}, these must be all > 0"
            )
            raise ValueError("Cannot resize the interpolation support.")
        self._origin = origin
        self.onGeometryChange()

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
        self._nsteps = np.ceil(length).astype(np.int64)
        self._nsteps[self._nsteps == 0] = 3
        if np.any(~(self._nsteps > 0)):
            logger.error(
                f"Cannot resize the interpolation support. The proposed number of steps is {self._nsteps}, these must be all > 0"
            )
            raise ValueError("Cannot resize the interpolation support.")
        self.onGeometryChange()

    @property
    def n_nodes(self):
        return np.prod(self.nsteps)

    @property
    def n_elements(self):
        return np.prod(self.nsteps_cells)

    @property
    def elements(self):
        global_index = np.arange(self.n_elements)
        cell_indexes = self.global_index_to_cell_index(global_index)

        return self.global_node_indices(self.cell_corner_indexes(cell_indexes))

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
        return np.array([xx.flatten(order="F"), yy.flatten(order="F"), zz.flatten(order="F")]).T

    def rotate(self, pos):
        """ """
        return np.einsum("ijk,ik->ij", self.rotation_xy[None, :, :], pos)

    def position_to_cell_index(self, pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Get the indexes (i,j,k) of a cell
        that a point is inside


        Parameters
        ----------
        pos : np.array
            Nx3 array of xyz locations

        Returns
        -------
        np.ndarray
            N,3 i,j,k indexes of the cell that the point is in
        """
        inside = self.inside(pos)
        pos = self.check_position(pos)
        cell_indexes = np.zeros((pos.shape[0], 3), dtype=int)

        x = pos[:, 0] - self.origin[None, 0]
        y = pos[:, 1] - self.origin[None, 1]
        z = pos[:, 2] - self.origin[None, 2]
        cell_indexes[inside, 0] = x[inside] // self.step_vector[None, 0]
        cell_indexes[inside, 1] = y[inside] // self.step_vector[None, 1]
        cell_indexes[inside, 2] = z[inside] // self.step_vector[None, 2]

        return cell_indexes, inside

    def position_to_cell_global_index(self, pos):
        ix, iy, iz = self.position_to_cell_index(pos)

    def inside(self, pos):
        # check whether point is inside box
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= pos[:, i] < self.maximum[None, i]
        return inside

    def check_position(self, pos: np.ndarray) -> np.ndarray:
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
        pos = np.array(pos)
        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            print("Position array needs to be a list of points or a point")
            raise ValueError("Position array needs to be a list of points or a point")
        return pos

    def _global_indicies(self, indexes: np.ndarray, nsteps: np.ndarray) -> np.ndarray:
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
        if indexes.shape[-1] != 3:
            raise ValueError("Last dimensions should be ijk indexing")
        original_shape = indexes.shape
        indexes = indexes.reshape(-1, 3)
        gi = (
            indexes[:, 0]
            + nsteps[None, 0] * indexes[:, 1]
            + nsteps[None, 0] * nsteps[None, 1] * indexes[:, 2]
        )
        return gi.reshape(original_shape[:-1])

    def cell_corner_indexes(self, cell_indexes: np.ndarray) -> np.ndarray:
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

        corner_indexes = np.zeros((cell_indexes.shape[0], 8, 3), dtype=int)

        xcorner = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1, 0, 0, 1, 1])
        zcorner = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        corner_indexes[:, :, 0] = (
            cell_indexes[:, None, 0] + corner_indexes[:, :, 0] + xcorner[None, :]
        )
        corner_indexes[:, :, 1] = (
            cell_indexes[:, None, 1] + corner_indexes[:, :, 1] + ycorner[None, :]
        )
        corner_indexes[:, :, 2] = (
            cell_indexes[:, None, 2] + corner_indexes[:, :, 2] + zcorner[None, :]
        )

        return corner_indexes

    def position_to_cell_corners(self, pos):
        """Get the global indices of the vertices (corners) of the cell containing each point.

        Parameters
        ----------
        pos : np.array
            (N, 3) array of xyz coordinates representing the positions of N points.
    
        Returns
        -------
        globalidx : np.array
            (N, 8) array of global indices corresponding to the 8 corner nodes of the cell 
            each point lies in. If a point lies outside the support, its corresponding entry 
            will be set to -1.
        inside : np.array
            (N,) boolean array indicating whether each point is inside the support domain.
        """
        cell_indexes, inside = self.position_to_cell_index(pos)
        corner_indexes = self.cell_corner_indexes(cell_indexes)
        globalidx = self.global_node_indices(corner_indexes)
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
        node_indexes = self.global_index_to_node_index(gi)
        return self.node_indexes_to_position(node_indexes), inside

    def node_indexes_to_position(self, node_indexes: np.ndarray) -> np.ndarray:
        original_shape = node_indexes.shape
        node_indexes = node_indexes.reshape((-1, 3))
        xyz = np.zeros((node_indexes.shape[0], 3), dtype=float)
        xyz[:, 0] = self.origin[0] + self.step_vector[0] * node_indexes[:, 0]
        xyz[:, 1] = self.origin[1] + self.step_vector[1] * node_indexes[:, 1]
        xyz[:, 2] = self.origin[2] + self.step_vector[2] * node_indexes[:, 2]
        xyz = xyz.reshape(original_shape)
        return xyz

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
        cell_indexes = np.zeros((global_index.shape[0], 3), dtype=int)
        cell_indexes[:, 0] = global_index % self.nsteps_cells[0, None]
        cell_indexes[:, 1] = global_index // self.nsteps_cells[0, None] % self.nsteps_cells[1, None]
        cell_indexes[:, 2] = (
            global_index // self.nsteps_cells[0, None] // self.nsteps_cells[1, None]
        )
        return cell_indexes

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
        original_shape = global_index.shape
        global_index = global_index.reshape((-1))
        local_indexes = np.zeros((global_index.shape[0], 3), dtype=int)
        local_indexes[:, 0] = global_index % self.nsteps[0, None]
        local_indexes[:, 1] = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        local_indexes[:, 2] = global_index // self.nsteps[0, None] // self.nsteps[1, None]
        return local_indexes.reshape(*original_shape, 3)

    def global_node_indices(self, indexes) -> np.ndarray:
        """
        Convert from node indexes to global node index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._global_indicies(indexes, self.nsteps)

    def global_cell_indices(self, indexes) -> np.ndarray:
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._global_indicies(indexes, self.nsteps_cells)

    @property
    def element_size(self):
        return np.prod(self.step_vector)

    @property
    def element_scale(self):
        # all elements are the same size
        return 1.0

    def vtk(self, node_properties={}, cell_properties={}):
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")

        from pyvista import CellType

        celltype = np.full(self.n_elements, CellType.VOXEL, dtype=np.uint8)
        elements = np.hstack(
            [np.zeros(self.elements.shape[0], dtype=int)[:, None] + 8, self.elements]
        )
        elements = elements.flatten()
        grid = pv.UnstructuredGrid(elements, celltype, self.nodes)
        for key, value in node_properties.items():
            grid[key] = value
        for key, value in cell_properties.items():
            grid.cell_arrays[key] = value
        return grid
