"""
Pure 2D regular grid geometry: origin/nsteps/step_vector indexing.
"""

import numpy as np
from typing import Tuple

from ..utils import getLogger

logger = getLogger(__name__)


class StructuredGrid2DGeometry:
    """A 2D regular grid defined by an origin, step vector and number of steps.

    Note: unlike :class:`StructuredGrid3DGeometry`, ``nsteps`` here is taken
    literally as a node count with no cells->nodes translation -- this matches
    the pre-existing behaviour of ``interpolators.supports.StructuredGrid2D``,
    which is intentionally not unified with the 3D convention.
    """

    dimension = 2

    def __init__(
        self,
        origin=None,
        nsteps=None,
        step_vector=None,
    ):
        """

        Parameters
        ----------
        origin - 2d list or numpy array
        nsteps - 2d list or numpy array of ints
        step_vector - 2d list or numpy array of int
        """
        if origin is None:
            origin = np.zeros(2)
        if nsteps is None:
            nsteps = np.array([10, 10])
        if step_vector is None:
            step_vector = np.ones(2)
        self.nsteps = np.ceil(np.array(nsteps)).astype(int)
        self.step_vector = np.array(step_vector)
        self.origin = np.array(origin)
        self.maximum = origin + self.nsteps * self.step_vector

        self.dim = 2
        self.nsteps_cells = self.nsteps - 1
        self.n_cell_x = self.nsteps[0] - 1
        self.n_cell_y = self.nsteps[1] - 1

    @property
    def nodes(self):
        max = self.origin + self.nsteps_cells * self.step_vector
        x = np.linspace(self.origin[0], max[0], self.nsteps[0])
        y = np.linspace(self.origin[1], max[1], self.nsteps[1])
        xx, yy = np.meshgrid(x, y, indexing="ij")
        return np.array([xx.flatten(order="F"), yy.flatten(order="F")]).T

    @property
    def n_nodes(self):
        return self.nsteps[0] * self.nsteps[1]

    @property
    def n_elements(self):
        return self.nsteps_cells[0] * self.nsteps_cells[1]

    @property
    def element_size(self):
        return np.prod(self.step_vector)

    @property
    def elements(self) -> np.ndarray:
        global_index = np.arange(self.n_elements)
        cell_indexes = self.global_index_to_cell_index(global_index)

        return self.global_node_indices(self.cell_corner_indexes(cell_indexes))

    def print_geometry(self):
        logger.info("Origin: %f %f %f" % (self.origin[0], self.origin[1], self.origin[2]))
        logger.info(
            "Cell size: %f %f %f" % (self.step_vector[0], self.step_vector[1], self.step_vector[2])
        )
        max = self.origin + self.nsteps_cells * self.step_vector
        logger.info("Max extent: %f %f %f" % (max[0], max[1], max[2]))

    def cell_centres(self, global_index: np.ndarray) -> np.ndarray:
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
        cell_indexes = self.global_index_to_cell_index(global_index)
        cell_centres = np.zeros((cell_indexes.shape[0], 2))

        cell_centres[:, 0] = (
            self.origin[None, 0]
            + self.step_vector[None, 0] * 0.5
            + self.step_vector[None, 0] * cell_indexes[:, 0]
        )
        cell_centres[:, 1] = (
            self.origin[None, 1]
            + self.step_vector[None, 1] * 0.5
            + self.step_vector[None, 1] * cell_indexes[:, 1]
        )
        return cell_centres

    def position_to_cell_index(self, pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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
        inside = self.inside(pos)
        cell_indexes = np.zeros((pos.shape[0], 2))
        cell_indexes[:, 0] = pos[:, 0] - self.origin[None, 0]
        cell_indexes[:, 1] = pos[:, 1] - self.origin[None, 1]
        cell_indexes /= self.step_vector[None, :]
        return cell_indexes.astype(int), inside

    def inside(self, pos: np.ndarray) -> np.ndarray:
        # check whether point is inside box
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(self.dim):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
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

        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            raise ValueError("Position array needs to be a list of points or a point")

        return pos

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
            gi = np.arange(self.n_nodes)
            indexes = self.global_index_to_node_index(gi)
            edge_mask = (
                (indexes[:, 0] > 0)
                & (indexes[:, 0] < self.nsteps[0] - 1)
                & (indexes[:, 1] > 0)
                & (indexes[:, 1] < self.nsteps[1] - 1)
            )
            indexes = indexes[edge_mask, :].T
        if indexes.ndim != 2:
            logger.error("indexes.ndim = %s, expected 2", indexes.ndim)
            return
        # determine which neighbours to return default is diagonals included.
        if mask is None:
            mask = np.array([[-1, 0, 1, -1, 0, 1, -1, 0, 1], [1, 1, 1, 0, 0, 0, -1, -1, -1]])
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return (neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1, :, :]).astype(
            np.int64
        )

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
        corner_indexes = np.zeros((cell_indexes.shape[0], 4, 2), dtype=np.int64)
        xcorner = np.array([0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1])
        corner_indexes[:, :, 0] = (
            cell_indexes[:, None, 0] + corner_indexes[:, :, 0] + xcorner[None, :]
        )
        corner_indexes[:, :, 1] = (
            cell_indexes[:, None, 1] + corner_indexes[:, :, 1] + ycorner[None, :]
        )
        return corner_indexes

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
        cell_indexes = np.zeros((global_index.shape[0], 2), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps_cells[0, None]
        cell_indexes[:, 1] = global_index // self.nsteps_cells[0, None] % self.nsteps_cells[1, None]
        return cell_indexes

    def global_index_to_node_index(self, global_index):
        cell_indexes = np.zeros((global_index.shape[0], 2), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps[0, None]
        cell_indexes[:, 1] = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        return cell_indexes

    def _global_indices(self, indexes: np.ndarray, nsteps: np.ndarray) -> np.ndarray:
        if len(indexes.shape) == 1:
            raise ValueError("Indexes must be a 2D array")
        if indexes.shape[-1] != 2:
            raise ValueError("Last dimension of cell indexing needs to be ijk indexing")
        original_shape = indexes.shape
        indexes = indexes.reshape(-1, 2)
        gi = indexes[:, 0] + nsteps[0] * indexes[:, 1]
        return gi.reshape(original_shape[:-1])

    def global_cell_indices(self, indexes: np.ndarray) -> np.ndarray:
        return self._global_indices(indexes, self.nsteps_cells)

    def global_node_indices(self, indexes: np.ndarray) -> np.ndarray:
        return self._global_indices(indexes, self.nsteps)

    def node_indexes_to_position(self, node_indexes: np.ndarray) -> np.ndarray:

        original_shape = node_indexes.shape
        node_indexes = node_indexes.reshape((-1, 2))
        xy = np.zeros((node_indexes.shape[0], 2), dtype=float)
        xy[:, 0] = self.origin[0] + self.step_vector[0] * node_indexes[:, 0]
        xy[:, 1] = self.origin[1] + self.step_vector[1] * node_indexes[:, 1]
        xy = xy.reshape(original_shape)
        return xy

    def position_to_cell_corners(self, pos):
        """Get the global indices of the vertices (corner) nodes of the cell containing each point.

        Parameters
        ----------
        pos : np.array
            (N, 2) array of xy coordinates representing the positions of N points.

        Returns
        -------
        globalidx : np.array
            (N, 4) array of global indices corresponding to the 4 corner nodes of the cell
            each point lies in. If a point lies outside the support, its corresponding entry
            will be set to -1.
        inside : np.array
            (N,) boolean array indicating whether each point is inside the support domain.
        """
        corner_index, inside = self.position_to_cell_index(pos)
        corners = self.cell_corner_indexes(corner_index)
        globalidx = self.global_node_indices(corners)
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

        node_indexes = self.global_index_to_node_index(gi.flatten())
        return self.node_indexes_to_position(node_indexes), inside

    def vtk(self, node_properties=None, cell_properties=None, z=0.0):
        """
        Create a vtk unstructured grid from the mesh
        """
        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}

        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for this functionality")

        from pyvista import CellType

        points = np.zeros((self.n_nodes, 3))
        points[:, :2] = self.nodes
        points[:, 2] = z
        celltype = np.full(self.n_elements, CellType.QUAD, dtype=np.uint8)
        vtk_elements = self.elements[:, [0, 1, 3, 2]]
        elements = np.hstack(
            [np.full((vtk_elements.shape[0], 1), 4, dtype=int), vtk_elements.astype(np.int64)]
        ).ravel()
        grid = pv.UnstructuredGrid(elements, celltype, points)
        for key, value in node_properties.items():
            grid.point_data[key] = value
        for key, value in cell_properties.items():
            grid.cell_data[key] = value
        return grid
