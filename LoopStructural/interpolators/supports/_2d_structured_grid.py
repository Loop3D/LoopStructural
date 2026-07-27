"""
Cartesian grid for fold interpolator

"""

import logging

import numpy as np
from . import SupportType
from ._base_support import BaseSupport
from LoopStructural.geometry import StructuredGrid2DGeometry
from typing import Dict, Tuple
from .._operator import Operator

logger = logging.getLogger(__name__)


class StructuredGrid2D(BaseSupport):
    """ """

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
        self.type = SupportType.StructuredGrid2D
        self._geom = StructuredGrid2DGeometry(origin=origin, nsteps=nsteps, step_vector=step_vector)
        self.properties = {}

        self.regions = {}
        self.regions["everywhere"] = np.ones(self.n_nodes).astype(bool)

    @property
    def origin(self):
        return self._geom.origin

    @property
    def nsteps(self):
        return self._geom.nsteps

    @property
    def nsteps_cells(self):
        return self._geom.nsteps_cells

    @property
    def step_vector(self):
        return self._geom.step_vector

    @property
    def maximum(self):
        return self._geom.maximum

    @property
    def dim(self):
        return self._geom.dim

    @property
    def n_cell_x(self):
        return self._geom.n_cell_x

    @property
    def n_cell_y(self):
        return self._geom.n_cell_y

    @property
    def nodes(self):
        return self._geom.nodes

    @property
    def n_nodes(self):
        return self._geom.n_nodes

    def set_nelements(self, nelements) -> int:
        raise NotImplementedError("Cannot set number of elements for 2D structured grid")

    @property
    def n_elements(self):
        return self._geom.n_elements

    @property
    def element_size(self):
        return self._geom.element_size

    @property
    def barycentre(self):
        return self.cell_centres(np.arange(self.n_elements))

    @property
    def elements(self) -> np.ndarray:
        return self._geom.elements

    def print_geometry(self):
        self._geom.print_geometry()

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
        return self._geom.cell_centres(global_index)

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
        return self._geom.position_to_cell_index(pos)

    def inside(self, pos: np.ndarray) -> np.ndarray:
        return self._geom.inside(pos)

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
        return self._geom.check_position(pos)

    def bilinear(self, local_coords: np.ndarray) -> np.ndarray:
        """
        returns the bilinear interpolation for the local coordinates
        Parameters
        ----------
        x - double, array of doubles
        y - double, array of doubles
        z - double, array of doubles

        Returns
        -------
        array of interpolation coefficients

        """

        return np.array(
            [
                (1 - local_coords[:, 0]) * (1 - local_coords[:, 1]),
                local_coords[:, 0] * (1 - local_coords[:, 1]),
                (1 - local_coords[:, 0]) * local_coords[:, 1],
                local_coords[:, 0] * local_coords[:, 1],
            ]
        ).T

    def position_to_local_coordinates(self, pos: np.ndarray) -> np.ndarray:
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
        local_coords = np.zeros(pos.shape)
        local_coords[:, 0] = (
            (pos[:, 0] - self.origin[None, 0]) % self.step_vector[None, 0]
        ) / self.step_vector[None, 0]
        local_coords[:, 1] = (
            (pos[:, 1] - self.origin[None, 1]) % self.step_vector[None, 1]
        ) / self.step_vector[None, 1]

        return local_coords

    def position_to_dof_coefs(self, pos: np.ndarray):
        """
        global posotion to interpolation coefficients
        Parameters
        ----------
        pos

        Returns
        -------

        """
        local_coords = self.position_to_local_coordinates(pos)
        weights = self.bilinear(local_coords)
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
        return self._geom.neighbour_global_indexes(mask=mask, **kwargs)

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
        return self._geom.cell_corner_indexes(cell_indexes)

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        return self._geom.global_index_to_cell_index(global_index)

    def global_index_to_node_index(self, global_index):
        return self._geom.global_index_to_node_index(global_index)

    def _global_indices(self, indexes: np.ndarray, nsteps: np.ndarray) -> np.ndarray:
        return self._geom._global_indices(indexes, nsteps)

    def global_cell_indices(self, indexes: np.ndarray) -> np.ndarray:
        return self._geom.global_cell_indices(indexes)

    def global_node_indices(self, indexes: np.ndarray) -> np.ndarray:
        return self._geom.global_node_indices(indexes)

    def node_indexes_to_position(self, node_indexes: np.ndarray) -> np.ndarray:
        return self._geom.node_indexes_to_position(node_indexes)

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
        return self._geom.position_to_cell_corners(pos)

    def evaluate_value(self, evaluation_points: np.ndarray, property_array: np.ndarray):
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

        v[inside, :] = self.position_to_dof_coefs(evaluation_points[inside, :])
        v[inside, :] *= property_array[idc[inside, :]]
        return np.sum(v, axis=1)

    def evaluate_gradient(self, evaluation_points, property_array):
        T = np.zeros((evaluation_points.shape[0], 2, 4))
        _vertices, T, elements, inside = self.get_element_gradient_for_location(evaluation_points)
        # indices = np.array([self.position_to_cell_index(evaluation_points)])
        # idc = self.global_indicies(indices.swapaxes(0,1))
        # print(idc)
        T[inside, 0, :] *= property_array[self.elements[elements[inside]]]
        T[inside, 1, :] *= property_array[self.elements[elements[inside]]]
        # T[inside, 2, :] *= self.properties[property_name][idc[inside, :]]
        return np.array([np.sum(T[:, 0, :], axis=1), np.sum(T[:, 1, :], axis=1)]).T

    def get_element_gradient_for_location(
        self, pos
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculates the gradient matrix at location pos
        :param pos: numpy array of location Nx3
        :return: Nx3x4 matrix
        """
        pos = np.asarray(pos)
        T = np.zeros((pos.shape[0], 2, 4))
        local_coords = self.position_to_local_coordinates(pos)
        vertices, inside = self.position_to_cell_corners(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements = self.global_cell_indices(elements)

        T[:, 0, 0] = -(1 - local_coords[:, 1])
        T[:, 0, 1] = 1 - local_coords[:, 1]
        T[:, 0, 2] = -local_coords[:, 1]
        T[:, 0, 3] = local_coords[:, 1]

        T[:, 1, 0] = -(1 - local_coords[:, 0])
        T[:, 1, 1] = -local_coords[:, 0]
        T[:, 1, 2] = 1 - local_coords[:, 0]
        T[:, 1, 3] = local_coords[:, 0]

        return vertices, T, elements, inside

    def get_element_for_location(
        self, pos: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

        vertices, inside = self.position_to_cell_vertices(pos)
        vertices = np.array(vertices)
        # print("ver", vertices.shape)
        # vertices = vertices.reshape((vertices.shape[1], 8, 3))
        elements, inside = self.position_to_cell_corners(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements = self.global_cell_indices(elements)
        a = self.position_to_dof_coefs(pos)
        return vertices, a, elements, inside

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
        return self._geom.position_to_cell_vertices(pos)

    def onGeometryChange(self):
        pass

    def vtk(self, node_properties=None, cell_properties=None, z=0.0):
        """
        Create a vtk unstructured grid from the mesh
        """
        return self._geom.vtk(node_properties=node_properties, cell_properties=cell_properties, z=z)

    def get_operators(self, weights: Dict[str, float]) -> Dict[str, Tuple[np.ndarray, float]]:
        """Get the finite difference mask operators used to build the smoothing/regularisation
        constraints for the 2d grid, scaled by the supplied weights.

        Parameters
        ----------
        weights : Dict[str, float]
            dictionary mapping operator name ("dxy", "dxx", "dyy") to its weighting factor

        Returns
        -------
        Dict[str, Tuple[np.ndarray, float]]
            dictionary mapping operator name to a tuple of (finite difference mask, weight)
        """
        # in a map we only want the xy operators
        operators = {
            "dxy": (Operator.Dxy_mask[1, :, :], weights["dxy"] * 2),
            "dxx": (Operator.Dxx_mask[1, :, :], weights["dxx"]),
            "dyy": (Operator.Dyy_mask[1, :, :], weights["dyy"]),
        }
        return operators
