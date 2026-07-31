"""Pure 3D regular grid geometry: origin/nsteps/step_vector indexing."""


import numpy as np

from loop_common.logging import get_logger as getLogger
from loop_common.utils import LoopException

logger = getLogger(__name__)


class StructuredGrid3DGeometry:
    """A 3D regular grid defined by an origin, step vector and number of steps."""

    dimension = 3

    def __init__(self, origin=None, nsteps=None, step_vector=None, rotation_xy=None):
        if origin is None:
            origin = np.zeros(3)
        if nsteps is None:
            nsteps = np.array([10, 10, 10])
        if step_vector is None:
            step_vector = np.ones(3)
        origin = np.array(origin)
        nsteps = np.array(nsteps)
        step_vector = np.array(step_vector)
        if np.any(step_vector == 0):
            logger.warning(f"Step vector {step_vector} has zero values")
        if np.any(nsteps == 0):
            raise LoopException("nsteps cannot be zero")
        if np.any(nsteps < 0):
            raise LoopException("nsteps cannot be negative")
        self._nsteps = np.array(nsteps, dtype=int)
        self._step_vector = np.array(step_vector)
        self._origin = np.array(origin)
        self._rotation_xy = np.zeros((3, 3))
        self._rotation_xy[0, 0] = 1
        self._rotation_xy[1, 1] = 1
        self._rotation_xy[2, 2] = 1
        self.rotation_xy = rotation_xy

    @property
    def volume(self):
        return np.prod(self.maximum - self.origin)

    def set_nelements(self, nelements) -> int:
        box_vol = self.volume
        ele_vol = box_vol / nelements
        step_vector = np.zeros(3)
        step_vector[:] = ele_vol ** (1.0 / 3.0)
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

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, nsteps):
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
        if rotation_xy is None:
            return
        if isinstance(rotation_xy, (float, int)):
            rotation_xy = np.array([[np.cos(np.deg2rad(rotation_xy)), -np.sin(np.deg2rad(rotation_xy)), 0], [np.sin(np.deg2rad(rotation_xy)), np.cos(np.deg2rad(rotation_xy)), 0], [0, 0, 1]])
        rotation_xy = np.array(rotation_xy)
        if rotation_xy.shape != (3, 3):
            raise ValueError(f"Rotation matrix should be 3x3, not {rotation_xy.shape}")
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

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        origin = np.array(origin)
        length = self.maximum - origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(np.int64)
        self._nsteps[self._nsteps == 0] = 3
        if np.any(~(self._nsteps > 0)):
            logger.error(f"Cannot resize the grid. The proposed number of steps is {self._nsteps}, these must be all > 0")
            raise ValueError("Cannot resize the grid.")
        self._origin = origin

    @property
    def maximum(self):
        return self.origin + self.nsteps_cells * self.step_vector

    @maximum.setter
    def maximum(self, maximum):
        maximum = np.array(maximum, dtype=float)
        length = maximum - self.origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(np.int64)
        self._nsteps[self._nsteps == 0] = 3
        if np.any(~(self._nsteps > 0)):
            logger.error(f"Cannot resize the grid. The proposed number of steps is {self._nsteps}, these must be all > 0")
            raise ValueError("Cannot resize the grid.")

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
            "LoopStructural grid geometry:  \n"
            f"Origin: {self.origin[0]} {self.origin[1]} {self.origin[2]} \n"
            f"Maximum: {self.maximum[0]} {self.maximum[1]} {self.maximum[2]} \n"
            f"Step Vector: {self.step_vector[0]} {self.step_vector[1]} {self.step_vector[2]} \n"
            f"Number of Steps: {self.nsteps[0]} {self.nsteps[1]} {self.nsteps[2]} \n"
            f"Degrees of freedon {self.n_nodes}"
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
        return np.einsum("ijk,ik->ij", self.rotation_xy[None, :, :], pos)

    def position_to_cell_index(self, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        inside = self.inside(pos)
        pos = self.check_position(pos)
        cell_indexes = np.zeros((pos.shape[0], 3), dtype=int)
        cell_indexes[:, 0] = (pos[:, 0] - self.origin[0]) / self.step_vector[0]
        cell_indexes[:, 1] = (pos[:, 1] - self.origin[1]) / self.step_vector[1]
        cell_indexes[:, 2] = (pos[:, 2] - self.origin[2]) / self.step_vector[2]
        return cell_indexes.astype(int), inside

    def inside(self, pos: np.ndarray) -> np.ndarray:
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[i]
            inside *= pos[:, i] < self.maximum[i]
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
                (indexes[:, 0] > 0) & (indexes[:, 0] < self.nsteps[0] - 1)
                & (indexes[:, 1] > 0) & (indexes[:, 1] < self.nsteps[1] - 1)
                & (indexes[:, 2] > 0) & (indexes[:, 2] < self.nsteps[2] - 1)
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
        corner_indexes = np.zeros((cell_indexes.shape[0], 8, 3), dtype=np.int64)
        xcorner = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1, 0, 0, 1, 1])
        zcorner = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        corner_indexes[:, :, 0] = cell_indexes[:, None, 0] + corner_indexes[:, :, 0] + xcorner[None, :]
        corner_indexes[:, :, 1] = cell_indexes[:, None, 1] + corner_indexes[:, :, 1] + ycorner[None, :]
        corner_indexes[:, :, 2] = cell_indexes[:, None, 2] + corner_indexes[:, :, 2] + zcorner[None, :]
        return corner_indexes

    def global_index_to_cell_index(self, global_index):
        cell_indexes = np.zeros((global_index.shape[0], 3), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps_cells[0, None]
        cell_indexes[:, 1] = (global_index // self.nsteps_cells[0, None]) % self.nsteps_cells[1, None]
        cell_indexes[:, 2] = (global_index // (self.nsteps_cells[0, None] * self.nsteps_cells[1, None])) % self.nsteps_cells[2, None]
        return cell_indexes

    def global_index_to_node_index(self, global_index):
        cell_indexes = np.zeros((global_index.shape[0], 3), dtype=np.int64)
        cell_indexes[:, 0] = global_index % self.nsteps[0, None]
        cell_indexes[:, 1] = (global_index // self.nsteps[0, None]) % self.nsteps[1, None]
        cell_indexes[:, 2] = (global_index // (self.nsteps[0, None] * self.nsteps[1, None])) % self.nsteps[2, None]
        return cell_indexes

    def global_node_indices(self, node_indexes):
        return node_indexes
