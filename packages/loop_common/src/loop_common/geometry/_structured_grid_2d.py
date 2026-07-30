"""Pure 2D regular grid geometry: origin/nsteps/step_vector indexing."""

import numpy as np
from typing import Tuple
from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


class StructuredGrid2DGeometry:
    """A 2D regular grid defined by an origin, step vector and number of steps."""

    dimension = 2

    def __init__(self, origin=None, nsteps=None, step_vector=None):
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
        inside = self.inside(pos)
        cell_indexes = np.zeros((pos.shape[0], 2))
        cell_indexes[:, 0] = pos[:, 0] - self.origin[None, 0]
        cell_indexes[:, 1] = pos[:, 1] - self.origin[None, 1]
        cell_indexes /= self.step_vector[None, :]
        return cell_indexes.astype(int), inside

    def inside(self, pos: np.ndarray) -> np.ndarray:
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(self.dim):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= pos[:, i] < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
        return inside

    def check_position(self, pos: np.ndarray) -> np.ndarray:
        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            raise ValueError("Position array needs to be a list of points or a point")
        return pos

    def neighbour_global_indexes(self, mask=None, **kwargs):
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
        if mask is None:
            mask = np.array([[-1, 0, 1, -1, 0, 1, -1, 0, 1], [1, 1, 1, 0, 0, 0, -1, -1, -1]])
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return (neighbours[0, :, :] + self.nsteps[0, None, None] * neighbours[1, :, :]).astype(np.int64)

    def cell_corner_indexes(self, cell_indexes: np.ndarray) -> np.ndarray:
        corner_indexes = np.zeros((cell_indexes.shape[0], 4, 2), dtype=np.int64)
        xcorner = np.array([0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1])
        corner_indexes[:, :, 0] = cell_indexes[:, None, 0] + corner_indexes[:, :, 0] + xcorner[None, :]
        corner_indexes[:, :, 1] = cell_indexes[:, None, 1] + corner_indexes[:, :, 1] + ycorner[None, :]
        return corner_indexes

    def global_index_to_cell_index(self, global_index):
        cell_indexes = np.zeros((global_index.shape[0], 2), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps_cells[0, None]
        cell_indexes[:, 1] = global_index // self.nsteps_cells[0, None] % self.nsteps_cells[1, None]
        return cell_indexes

    def global_index_to_node_index(self, global_index):
        cell_indexes = np.zeros((global_index.shape[0], 2), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps[0, None]
        cell_indexes[:, 1] = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        return cell_indexes

    def global_node_indices(self, node_indexes):
        return node_indexes
