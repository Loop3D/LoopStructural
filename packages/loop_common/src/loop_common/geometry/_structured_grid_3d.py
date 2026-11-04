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

    def _on_geometry_change(self):
        """Hook called after origin, step_vector, nsteps or maximum change; subclasses override."""

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
        self._on_geometry_change()

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
        self._on_geometry_change()

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
        self._on_geometry_change()

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
        self._on_geometry_change()

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
        """Get the (i, j, k) index of the cell that each point is in.

        Points are expected in the grid's own frame; ``rotation_xy`` is not applied
        (see :meth:`rotate`).

        Parameters
        ----------
        pos : np.ndarray
            Nx3 array of xyz locations

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            N,3 integer cell indexes and a boolean mask of the points that are inside the
            grid. Points on the maximum boundary belong to the last cell. The index of a point
            outside the grid is (0, 0, 0): use the mask, not the index, to detect them.
        """
        pos = self.check_position(pos)
        inside = self.inside(pos)
        cell_indexes = np.zeros((pos.shape[0], 3), dtype=int)
        for axis in range(3):
            cell = np.floor((pos[inside, axis] - self.origin[axis]) / self.step_vector[axis])
            cell_indexes[inside, axis] = np.clip(cell.astype(int), 0, self.nsteps_cells[axis] - 1)
        return cell_indexes, inside

    def inside(self, pos: np.ndarray) -> np.ndarray:
        """Boolean mask of the points inside the grid, boundary included (to a small tolerance)."""
        pos = self.check_position(pos)
        extent = np.maximum(np.abs(self.maximum - self.origin), 1.0)
        tol = np.finfo(float).eps * 32.0 * extent
        return np.all((pos >= (self.origin - tol)) & (pos <= (self.maximum + tol)), axis=1)

    def check_position(self, pos: np.ndarray) -> np.ndarray:
        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            raise ValueError("Position array needs to be a list of points or a point")
        return pos

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
        cell_indexes[:, 2] = global_index // self.nsteps_cells[0, None] // self.nsteps_cells[1, None]
        return cell_indexes

    def global_index_to_node_index(self, global_index):
        """Node ijk indexes (shape ``(*global_index.shape, 3)``) from flat node indexes."""
        global_index = np.asarray(global_index)
        original_shape = global_index.shape
        global_index = global_index.reshape(-1)
        node_indexes = np.zeros((global_index.shape[0], 3), dtype=np.int64)
        node_indexes[:, 0] = global_index % self.nsteps[0, None]
        node_indexes[:, 1] = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        node_indexes[:, 2] = global_index // self.nsteps[0, None] // self.nsteps[1, None]
        return node_indexes.reshape(*original_shape, 3)

    def _global_indicies(self, indexes: np.ndarray, nsteps: np.ndarray) -> np.ndarray:
        """Convert (..., 3) ijk indexes to flat indexes for a grid with ``nsteps`` per axis (x fastest)."""
        indexes = np.asarray(indexes)
        if indexes.ndim == 1:
            raise ValueError("Cell indexes needs to be Nx3")
        if indexes.shape[-1] != 3:
            raise ValueError("Last dimensions should be ijk indexing")
        original_shape = indexes.shape
        indexes = indexes.reshape(-1, 3)
        gi = indexes[:, 0] + nsteps[0] * indexes[:, 1] + nsteps[0] * nsteps[1] * indexes[:, 2]
        return gi.reshape(original_shape[:-1])

    def global_node_indices(self, node_indexes):
        """Flat node indexes from (..., 3) node ijk indexes."""
        return self._global_indicies(node_indexes, self.nsteps)

    def global_cell_indices(self, cell_indexes):
        """Flat cell indexes from (..., 3) cell ijk indexes; the order cell data is stored in."""
        return self._global_indicies(cell_indexes, self.nsteps_cells)

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
        nx, ny = self.nsteps[0], self.nsteps[1]
        offsets = np.array(
            [0, 1, nx, nx + 1, nx * ny, nx * ny + 1, nx * ny + nx, nx * ny + nx + 1],
            dtype=np.intp,
        )
        g = cell_indexes[:, 0] + nx * cell_indexes[:, 1] + nx * ny * cell_indexes[:, 2]
        globalidx = g[:, None] + offsets[None, :]  # (N, 8)
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

    @property
    def element_size(self):
        return np.prod(self.step_vector)

    @property
    def element_scale(self):
        # all elements are the same size
        return 1.0

    def trilinear(self, local_coords):
        """
        Returns the trilinear interpolation for the local coordinates

        Parameters
        ----------
        x - double, array of doubles
        y - double, array of doubles
        z - double, array of doubles

        Returns
        -------
        array of interpolation coefficients

        """
        interpolant = np.zeros((local_coords.shape[0], 8), dtype=np.float64)
        # interpolant[:,0] = (1 - local_coords[:, 0])* (1 - local_coords[:, 1])* (1 - local_coords[:, 2])
        # interpolant[:,1] = local_coords[:, 0]* (1 - local_coords[:, 1])* (1 - local_coords[:, 2])
        # interpolant[:,2] = (1 - local_coords[:, 0])* local_coords[:, 1]* (1 - local_coords[:, 2])
        # interpolant[:,3] = (1 - local_coords[:, 0])* (1 - local_coords[:, 1])* local_coords[:, 2]
        # interpolant[:,4] = local_coords[:, 0] * (1 - local_coords[:, 1]) * local_coords[:, 2]
        # interpolant[:,5] = (1 - local_coords[:, 0]) * local_coords[:, 1] * local_coords[:, 2]
        # interpolant[:,6] = local_coords[:, 0] * local_coords[:, 1] * (1 - local_coords[:, 2])
        # interpolant[:,7] =  local_coords[:, 0] * local_coords[:, 1] * local_coords[:, 2]
        interpolant[:, 0] = (
            (1 - local_coords[:, 0]) * (1 - local_coords[:, 1]) * (1 - local_coords[:, 2])
        )
        interpolant[:, 1] = local_coords[:, 0] * (1 - local_coords[:, 1]) * (1 - local_coords[:, 2])
        interpolant[:, 2] = (1 - local_coords[:, 0]) * local_coords[:, 1] * (1 - local_coords[:, 2])
        interpolant[:, 4] = (1 - local_coords[:, 0]) * (1 - local_coords[:, 1]) * local_coords[:, 2]
        interpolant[:, 5] = local_coords[:, 0] * (1 - local_coords[:, 1]) * local_coords[:, 2]
        interpolant[:, 6] = (1 - local_coords[:, 0]) * local_coords[:, 1] * local_coords[:, 2]
        interpolant[:, 3] = local_coords[:, 0] * local_coords[:, 1] * (1 - local_coords[:, 2])
        interpolant[:, 7] = local_coords[:, 0] * local_coords[:, 1] * local_coords[:, 2]
        return interpolant

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
        pos = self.check_position(pos)
        cell_indexes, _inside = self.position_to_cell_index(pos)

        # Use cell origin (not modulo) so points on max faces/vertices map to local=1.
        cell_origin = self.origin[None, :] + self.step_vector[None, :] * cell_indexes
        local_coords = (pos - cell_origin) / self.step_vector[None, :]
        local_coords = np.clip(local_coords, 0.0, 1.0)
        return local_coords

    def position_to_dof_coefs(self, pos):
        """
        Global posotion to interpolation coefficients

        Parameters
        ----------
        pos

        Returns
        -------

        """
        local_coords = self.position_to_local_coordinates(pos)
        weights = self.trilinear(local_coords)
        return weights

    def evaluate_value(self, evaluation_points, property_array):
        """Trilinear interpolation of node values at the evaluation points.

        Parameters
        ----------
        evaluation_points : np.ndarray
            N x 3 locations in the grid's own frame
        property_array : np.ndarray
            values at the nodes, shape ``(n_nodes,)`` or ``(n_nodes, k)``, in node order
            (x varies fastest)

        Returns
        -------
        np.ndarray
            ``(N,)`` (or ``(N, k)``) interpolated values; NaN for points outside the grid
        """
        property_array = np.asarray(property_array)
        evaluation_points = self.check_position(evaluation_points)
        if property_array.shape[0] != self.n_nodes:
            logger.error("Property array does not match grid")
            raise ValueError(
                f"cannot assign {property_array.shape[0]} vlaues to array of shape {self.n_nodes}"
            )
        idc, inside = self.position_to_cell_corners(evaluation_points)
        if idc.shape[0] != inside.shape[0]:
            raise ValueError("index does not match number of nodes")
        weights = self.position_to_dof_coefs(evaluation_points[inside, :])
        if property_array.ndim == 1:
            v = np.zeros(idc.shape)
            v[:, :] = np.nan
            v[inside, :] = weights
            v[inside, :] *= property_array[idc[inside, :]]
            return np.sum(v, axis=1)
        result = np.full((idc.shape[0],) + property_array.shape[1:], np.nan)
        result[inside] = np.einsum("nc,nc...->n...", weights, property_array[idc[inside, :]])
        return result
