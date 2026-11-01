"""
Cartesian grid for fold interpolator

"""

import logging

import numpy as np

from loop_common.geometry import StructuredGrid2DGeometry

from ..math.finite_difference_stencil import Operator
from . import SupportType
from ._base_support import BaseSupport

logger = logging.getLogger(__name__)


class StructuredGrid2D(StructuredGrid2DGeometry, BaseSupport):
    """A 2D structured grid used as an interpolation support.

    The grid geometry (node and cell positions, point location, index conversions, bilinear
    weights and node-value interpolation) comes from
    :class:`loop_common.geometry.StructuredGrid2DGeometry`; this class adds the interpolation
    machinery. ``nsteps`` is the number of **nodes** in each direction.
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
        StructuredGrid2DGeometry.__init__(self, origin=origin, nsteps=nsteps, step_vector=step_vector)
        self.type = SupportType.StructuredGrid2D
        self.properties = {}
        self.regions = {}
        self.regions["everywhere"] = np.ones(self.n_nodes).astype(bool)

    def set_nelements(self, nelements) -> int:
        raise NotImplementedError("Cannot set number of elements for 2D structured grid")

    @property
    def barycentre(self):
        return self.cell_centres(np.arange(self.n_elements))

    # @property
    # def barycentre(self):
    #     return self.cell_centres(np.arange(self.n_elements))
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
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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

        T[:, 0, :] /= self.step_vector[None, 0]
        T[:, 1, :] /= self.step_vector[None, 1]
        return vertices, T, elements, inside

    def get_element_for_location(
        self, pos: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

        vertices, inside = self.position_to_cell_vertices(pos)
        vertices = np.array(vertices)
        # print("ver", vertices.shape)
        # vertices = vertices.reshape((vertices.shape[1], 8, 3))
        elements, inside = self.position_to_cell_corners(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements = self.global_cell_indices(elements)
        a = self.position_to_dof_coefs(pos)
        return vertices, a, elements, inside

    def onGeometryChange(self):
        pass

    def vtk(self, z, *, node_properties=None, cell_properties=None):
        """
        Create a vtk unstructured grid from the mesh
        """
        import pyvista as pv
        from pyvista import CellType

        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
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

    def get_operators(self, weights: dict[str, float]) -> dict[str, tuple[np.ndarray, float]]:
        """Get

        Parameters
        ----------
        weights : Dict[str, float]
            _description_

        Returns
        -------
        Dict[str, Tuple[np.ndarray, float]]
            _description_
        """
        # in a map we only want the xy operators
        operators = {
            "dxy": (Operator.Dxy_mask[1, :, :], weights["dxy"] * 2),
            "dxx": (Operator.Dxx_mask[1, :, :], weights["dxx"]),
            "dyy": (Operator.Dyy_mask[1, :, :], weights["dyy"]),
        }
        return operators
