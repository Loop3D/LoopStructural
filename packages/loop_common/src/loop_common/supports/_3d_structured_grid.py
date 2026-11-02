"""
Cartesian grid for fold interpolator

"""

import numpy as np

from loop_common.logging import get_logger as getLogger

from ..math.finite_difference_stencil import Operator
from . import SupportType
from ._3d_base_structured import BaseStructuredSupport

logger = getLogger(__name__)


class StructuredGrid(BaseStructuredSupport):
    """ """

    def __init__(
        self,
        origin=None,
        nsteps_cells=None,
        step_vector=None,
        nsteps=None,
        rotation_xy=None,
        properties=None,
        cell_properties=None,
        name="StructuredGrid",
    ):
        """

        Parameters
        ----------
        origin : array-like, optional
            Grid origin.
        nsteps : array-like, optional
            Legacy node-count alias. Treated as cell counts for compatibility.
        nsteps_cells : array-like, optional
            Number of cells in each direction.
        step_vector : array-like, optional
            Cell size in each direction.
        """
        if origin is None:
            origin = np.zeros(3)
        if step_vector is None:
            step_vector = np.ones(3)
        if nsteps_cells is None:
            nsteps_cells = np.array([10, 10, 10]) if nsteps is None else np.array(nsteps)
        BaseStructuredSupport.__init__(
            self,
            origin,
            nsteps_cells,
            step_vector,
            rotation_xy=rotation_xy,
        )
        self.type = SupportType.StructuredGrid
        self.regions = {}
        self.regions["everywhere"] = np.ones(self.n_nodes).astype(bool)

    def onGeometryChange(self):
        if self.interpolator is not None:
            self.interpolator.reset()

    @property
    def barycentre(self):
        return self.cell_centres(np.arange(self.n_elements))

    def cell_centres(self, global_index):
        """Get the centre of specified cells

        Parameters
        ----------
        global_index : array/list
            container of integer global indexes to cells

        Returns
        -------
        numpy array
            Nx3 array of cell centres
        """
        cell_indexes = self.global_index_to_cell_index(global_index)
        return (
            self.origin[None, :]
            + self.step_vector[None, :] * (cell_indexes)
            + self.step_vector[None, :] * 0.5
        )

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
                & (indexes[:, 2] > 0)
                & (indexes[:, 2] < self.nsteps[2] - 1)
            )
            indexes = indexes[edge_mask, :].T

        # indexes = np.array(indexes).T
        if indexes.ndim != 2:
            return
        # determine which neighbours to return default is diagonals included.
        if mask is None:
            mask = np.array(
                [
                    [
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                        -1,
                        0,
                        1,
                    ],
                    [
                        -1,
                        -1,
                        -1,
                        0,
                        0,
                        0,
                        1,
                        1,
                        1,
                        -1,
                        -1,
                        -1,
                        0,
                        0,
                        0,
                        1,
                        1,
                        1,
                        -1,
                        -1,
                        -1,
                        0,
                        0,
                        0,
                        1,
                        1,
                        1,
                    ],
                    [
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                    ],
                ]
            )
        neighbours = indexes[:, None, :] + mask[:, :, None]
        return (
            neighbours[0, :, :]
            + self.nsteps[0, None, None] * neighbours[1, :, :]
            + self.nsteps[0, None, None] * self.nsteps[1, None, None] * neighbours[2, :, :]
        ).astype(np.int64)

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
            raise ValueError(
                f"cannot assign {property_array.shape[0]} vlaues to array of shape {self.n_nodes}"
            )
        if self._has_uniform_lookup():
            idc, inside, local_coords = self._lookup_points(evaluation_points)
            coefs = self.trilinear(local_coords[inside, :])
        else:
            idc, inside = self.position_to_cell_corners(evaluation_points)
            coefs = self.position_to_dof_coefs(evaluation_points[inside, :])
        v = np.zeros(idc.shape)
        v[:, :] = np.nan
        v[inside, :] = coefs
        v[inside, :] *= property_array[idc[inside, :]]
        return np.sum(v, axis=1)

    def evaluate_gradient(self, evaluation_points, property_array) -> np.ndarray:
        """Evaluate the gradient at a location given node values

        Parameters
        ----------
        evaluation_points : np.array((N,3))
            locations
        property_array : np.array((self.nx))
            value node, has to be the same length as the number of nodes

        Returns
        -------
        np.array((N,3),dtype=float)
            gradient of the implicit function at the locations

        Raises
        ------
        ValueError
            if the array is not the same shape as the number of nodes

        Notes
        -----
        The implicit function gradient is not normalised, to convert to
        a unit vector normalise using vector/=np.linalg.norm(vector,axis=1)[:,None]
        """
        if property_array.shape[0] != self.n_nodes:
            logger.error("Property array does not match grid")
            raise ValueError(
                f"cannot assign {property_array.shape[0]} vlaues to array of shape {self.n_nodes}"
            )

        if self._has_uniform_lookup():
            evaluation_points = self.check_position(evaluation_points)
            idc, inside, local_coords = self._lookup_points(evaluation_points)
            weights = self._gradient_weights(local_coords[inside, :])
        else:
            idc, inside = self.position_to_cell_corners(evaluation_points)
            weights = self.get_element_gradient_for_location(evaluation_points[inside, :])[1]
        T = np.zeros((idc.shape[0], 3, 8))
        T[inside, :, :] = weights
        if np.any(inside) and np.max(idc[inside, :]) > property_array.shape[0]:
            cix, ciy, ciz = self.position_to_cell_index(evaluation_points)
            if not np.all(cix[inside] < self.nsteps_cells[0]):
                print(
                    evaluation_points[inside, :][cix[inside] < self.nsteps_cells[0], 0],
                    self.origin[0],
                    self.maximum[0],
                )
            if not np.all(ciy[inside] < self.nsteps_cells[1]):
                print(
                    evaluation_points[inside, :][ciy[inside] < self.nsteps_cells[1], 1],
                    self.origin[1],
                    self.maximum[1],
                )
            if not np.all(ciz[inside] < self.nsteps_cells[2]):
                print(ciz[inside], self.nsteps_cells[2])
                print(self.step_vector, self.nsteps_cells, self.nsteps)
                print(
                    evaluation_points[inside, :][~(ciz[inside] < self.nsteps_cells[2]), 2],
                    self.origin[2],
                    self.maximum[2],
                )

            raise ValueError("index does not match number of nodes")
        T[inside, 0, :] *= property_array[idc[inside, :]]
        T[inside, 1, :] *= property_array[idc[inside, :]]
        T[inside, 2, :] *= property_array[idc[inside, :]]
        return np.array(
            [
                np.sum(T[:, 0, :], axis=1),
                np.sum(T[:, 1, :], axis=1),
                np.sum(T[:, 2, :], axis=1),
            ]
        ).T

    def _has_uniform_lookup(self):
        """Check that the cell geometry is the uniform geometry of this class.

        A subclass with its own local coordinates (for example a grid with
        different cell sizes) cannot use _lookup_points and _gradient_weights.
        """
        cls = type(self)
        return (
            cls.position_to_local_coordinates is StructuredGrid.position_to_local_coordinates
            and cls.get_element_gradient_for_location
            is StructuredGrid.get_element_gradient_for_location
            and cls.position_to_cell_corners is StructuredGrid.position_to_cell_corners
        )

    def _lookup_points(self, pos):
        """Find the cell, the corner nodes and the local coordinates of points.

        The cell index is found one time for all points.

        Parameters
        ----------
        pos : np.array((N,3),dtype=float)
            locations

        Returns
        -------
        corners, inside, local_coords
            corners is an integer (N,8) array of node indexes, with -1 in each
            row of a point that is outside. inside is a boolean (N,) array.
            local_coords is a float (N,3) array in the range 0 to 1. Rows of
            points that are outside use the cell with index 0.
        """
        pos = self.check_position(pos)
        cell_indexes, inside = self.position_to_cell_index(pos)
        nx, ny = self.nsteps[0], self.nsteps[1]
        offsets = np.array(
            [0, 1, nx, nx + 1, nx * ny, nx * ny + 1, nx * ny + nx, nx * ny + nx + 1],
            dtype=np.intp,
        )
        g = cell_indexes[:, 0] + nx * cell_indexes[:, 1] + nx * ny * cell_indexes[:, 2]
        corners = g[:, None] + offsets[None, :]
        corners[~inside] = -1
        # Use cell origin (not modulo) so points on max faces/vertices map to local=1.
        cell_origin = self.origin[None, :] + self.step_vector[None, :] * cell_indexes
        local_coords = (pos - cell_origin) / self.step_vector[None, :]
        local_coords = np.clip(local_coords, 0.0, 1.0)
        return corners, inside, local_coords

    def _gradient_weights(self, local_coords):
        """Get the gradient weights of the 8 corner nodes of a cell.

        Parameters
        ----------
        local_coords : np.array((N,3),dtype=float)
            local coordinates of the points in their cells

        Returns
        -------
        np.array((N,3,8),dtype=float)
            weights for each axis and corner, divided by the step vector
        """
        T = np.zeros((local_coords.shape[0], 3, 8))
        T[:, 0, 0] = (1 - local_coords[:, 2]) * (local_coords[:, 1] - 1)  # v000
        T[:, 0, 1] = (1 - local_coords[:, 1]) * (1 - local_coords[:, 2])
        T[:, 0, 2] = -local_coords[:, 1] * (1 - local_coords[:, 2])
        T[:, 0, 4] = -(1 - local_coords[:, 1]) * local_coords[:, 2]
        T[:, 0, 5] = (1 - local_coords[:, 1]) * local_coords[:, 2]
        T[:, 0, 6] = -local_coords[:, 1] * local_coords[:, 2]
        T[:, 0, 3] = local_coords[:, 1] * (1 - local_coords[:, 2])
        T[:, 0, 7] = local_coords[:, 1] * local_coords[:, 2]

        T[:, 1, 0] = (local_coords[:, 0] - 1) * (1 - local_coords[:, 2])
        T[:, 1, 1] = -local_coords[:, 0] * (1 - local_coords[:, 2])
        T[:, 1, 2] = (1 - local_coords[:, 0]) * (1 - local_coords[:, 2])
        T[:, 1, 4] = -(1 - local_coords[:, 0]) * local_coords[:, 2]
        T[:, 1, 5] = -local_coords[:, 0] * local_coords[:, 2]
        T[:, 1, 6] = (1 - local_coords[:, 0]) * local_coords[:, 2]
        T[:, 1, 3] = local_coords[:, 0] * (1 - local_coords[:, 2])
        T[:, 1, 7] = local_coords[:, 0] * local_coords[:, 2]

        T[:, 2, 0] = -(1 - local_coords[:, 0]) * (1 - local_coords[:, 1])
        T[:, 2, 1] = -local_coords[:, 0] * (1 - local_coords[:, 1])
        T[:, 2, 2] = -(1 - local_coords[:, 0]) * local_coords[:, 1]
        T[:, 2, 4] = (1 - local_coords[:, 0]) * (1 - local_coords[:, 1])
        T[:, 2, 5] = local_coords[:, 0] * (1 - local_coords[:, 1])
        T[:, 2, 6] = (1 - local_coords[:, 0]) * local_coords[:, 1]
        T[:, 2, 3] = -local_coords[:, 0] * local_coords[:, 1]
        T[:, 2, 7] = local_coords[:, 0] * local_coords[:, 1]
        T[:, 0, :] /= self.step_vector[None, 0]
        T[:, 1, :] /= self.step_vector[None, 1]
        T[:, 2, :] /= self.step_vector[None, 2]
        return T

    def get_element_gradient_for_location(self, pos: np.ndarray):
        """
        Get the gradient of the element at the locations.

        Parameters
        ----------
        pos : np.array((N,3),dtype=float)
            locations

        Returns
        -------
        vertices, gradient, element, inside
            [description]
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
        pos = np.asarray(pos)
        local_coords = self.position_to_local_coordinates(pos)
        vertices, inside = self.position_to_cell_vertices(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements = self.global_cell_indices(elements)

        T = self._gradient_weights(local_coords)
        return vertices, T, elements, inside

    def get_element_for_location(self, pos: np.ndarray):
        """Calculate the shape function of elements
        for a location

        Parameters
        ----------
        pos : np.array((N,3))
            location of points to calculate the shape function

        Returns
        -------
        [type]
            [description]
        """
        vertices, inside = self.position_to_cell_vertices(pos)
        vertices = np.array(vertices)
        # print("ver", vertices.shape)
        # vertices = vertices.reshape((vertices.shape[1], 8, 3))
        elements, inside = self.position_to_cell_corners(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements = self.global_cell_indices(elements)
        a = self.position_to_dof_coefs(pos)
        return vertices, a, elements, inside

    def get_elements(self):
        return

    def to_dict(self):
        return {
            "type": self.type.numerator,
            **super().to_dict(),
        }

    def get_operators(self, weights: dict[str, float]) -> dict[str, tuple[np.ndarray, float]]:
        """Gets the operators specific to this support

        Parameters
        ----------
        weights : Dict[str, float]
            weight value per operator

        Returns
        -------
        operators
            A dictionary with a numpy array and float weight
        """
        operators = {
            "dxy": (Operator.Dxy_mask, weights["dxy"] / 4),
            "dyz": (Operator.Dyz_mask, weights["dyz"] / 4),
            "dxz": (Operator.Dxz_mask, weights["dxz"] / 4),
            "dxx": (Operator.Dxx_mask, weights["dxx"] / 1),
            "dyy": (Operator.Dyy_mask, weights["dyy"] / 1),
            "dzz": (Operator.Dzz_mask, weights["dzz"] / 1),
        }
        return operators
