"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

from abc import abstractmethod
import logging
from typing import Tuple
import numpy as np
from scipy import sparse

from . import SupportType
from ._2d_structured_grid import StructuredGrid2D
from ._base_support import BaseSupport
from ._aabb import _initialise_aabb
from ._face_table import _init_face_table

logger = logging.getLogger(__name__)


class BaseUnstructured2d(BaseSupport):
    """ """

    dimension = 2

    def __init__(self, elements, vertices, neighbours, aabb_nsteps=None):
        self.type = SupportType.BaseUnstructured2d
        self._elements = elements
        self.vertices = vertices
        if self.elements.shape[1] == 3:
            self.order = 1
        elif self.elements.shape[1] == 6:
            self.order = 2
        self.dof = self.vertices.shape[0]
        self.neighbours = neighbours
        self.minimum = np.min(self.nodes, axis=0)
        self.maximum = np.max(self.nodes, axis=0)
        length = self.maximum - self.minimum
        self.minimum -= length * 0.1
        self.maximum += length * 0.1
        if aabb_nsteps is None:
            box_vol = np.prod(self.maximum - self.minimum)
            element_volume = box_vol / (len(self.elements) / 20)
            # calculate the step vector of a regular cube
            step_vector = np.zeros(2)
            step_vector[:] = element_volume ** (1.0 / 2.0)
            # number of steps is the length of the box / step vector
            aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(int)
            # make sure there is at least one cell in every dimension
            aabb_nsteps[aabb_nsteps < 2] = 2
        step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
        self.aabb_grid = StructuredGrid2D(self.minimum, nsteps=aabb_nsteps, step_vector=step_vector)
        # make a big table to store which tetra are in which element.
        # if this takes up too much memory it could be simplified by using sparse matrices or dict but
        # at the expense of speed
        self._aabb_table = sparse.csr_matrix(
            (self.aabb_grid.n_elements, len(self.elements)), dtype=bool
        )
        self._shared_element_relationships = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int
        )
        self._shared_elements = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], self.dimension), dtype=int
        )

    @property
    def aabb_table(self):
        if np.sum(self._aabb_table) == 0:
            _initialise_aabb(self)
        return self._aabb_table

    def set_nelements(self, nelements) -> int:
        raise NotImplementedError

    @property
    def shared_elements(self):
        if np.sum(self._shared_elements) == 0:
            _init_face_table(self)
        return self._shared_elements

    @property
    def shared_element_relationships(self):
        if np.sum(self._shared_element_relationships) == 0:
            _init_face_table(self)
        return self._shared_element_relationships

    @property
    def elements(self):
        return self._elements

    def onGeometryChange(self):
        pass

    @property
    def n_elements(self):
        return self.elements.shape[0]

    @property
    def n_nodes(self):
        return self.vertices.shape[0]

    def inside(self, pos):
        if pos.shape[1] > self.dimension:
            logger.warning(f"Converting {pos.shape[1]} to 3d using first {self.dimension} columns")
            pos = pos[:, : self.dimension]

        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(self.dimension):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= pos[:, i] < self.maximum[None, i]
        return inside

    @property
    def ncps(self):
        """
        Returns the number of nodes for an element in the mesh
        """
        return self.elements.shape[1]

    @property
    def nodes(self):
        """
        Gets the nodes of the mesh as a property rather than using a function, accessible as a property! Python magic!

        Returns
        -------
        nodes : np.array((N,3))
            Fortran ordered
        """
        return self.vertices

    @property
    def barycentre(self):
        """
        Return the barycentres of all tetrahedrons or of specified tetras using
        global index

        Parameters
        ----------
        elements - numpy array
            global index

        Returns
        -------

        """
        element_idx = np.arange(0, self.n_elements)
        elements = self.elements[element_idx]
        barycentre = np.sum(self.nodes[elements][:, :3, :], axis=1) / 3.0
        return barycentre

    @property
    def shared_element_norm(self):
        """
        Get the normal to all of the shared elements
        """
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        norm = np.zeros_like(v1)
        norm[:, 0] = v1[:, 1]
        norm[:, 1] = -v1[:, 0]
        return norm

    @property
    def shared_element_size(self):
        """
        Get the size of the shared elements
        """
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        return np.linalg.norm(v1, axis=1)

    @property
    def element_size(self):
        v1 = self.nodes[self.elements[:, 1], :] - self.nodes[self.elements[:, 0], :]
        v2 = self.nodes[self.elements[:, 2], :] - self.nodes[self.elements[:, 0], :]
        # cross product isn't defined in 2d, numpy returns the magnitude of the orthogonal vector.
        return 0.5 * np.cross(v1, v2, axisa=1, axisb=1)

    @abstractmethod
    def evaluate_shape(self, locations) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate the shape functions at the locations

        Parameters
        ----------
        locations - numpy array
            locations to evaluate

        Returns
        -------

        """
        pass

    def element_area(self, elements):
        tri_points = self.nodes[self.elements[elements, :], :]
        M_t = np.ones((tri_points.shape[0], 3, 3))
        M_t[:, :, 1:] = tri_points[:, :3, :]
        area = np.abs(np.linalg.det(M_t)) * 0.5
        return area

    def evaluate_value(self, evaluation_points: np.ndarray, property_array: np.ndarray):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - numpy array
            property values at nodes

        Returns
        -------

        """
        pos = np.asarray(evaluation_points)
        return_values = np.zeros(pos.shape[0])
        return_values[:] = np.nan
        _verts, c, tri, inside = self.get_element_for_location(pos[:, :2])
        inside = tri >= 0
        # vertices, c, elements, inside = self.get_elements_for_location(pos)
        return_values[inside] = np.sum(
            c[inside, :] * property_array[self.elements[tri[inside], :]], axis=1
        )
        return return_values

    def evaluate_gradient(self, evaluation_points, property_array):
        """
        Evaluate the gradient of an interpolant at the locations

        Parameters
        ----------
        pos - numpy array
            locations
        prop - string
            property to evaluate


        Returns
        -------

        """
        values = np.zeros(evaluation_points.shape)
        values[:] = np.nan
        element_gradients, tri, inside = self.evaluate_shape_derivatives(evaluation_points[:, :2])
        inside = tri >= 0

        values[inside, :] = (
            element_gradients[inside, :, :] * property_array[self.elements[tri[inside], :, None]]
        ).sum(1)
        return values

    def get_element_for_location(
        self,
        points: np.ndarray,
        return_verts=True,
        return_bc=True,
        return_inside=True,
        return_tri=True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Determine the elements from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        if return_verts:
            verts = np.zeros((points.shape[0], self.dimension + 1, self.dimension))
        else:
            verts = np.zeros((0, 0, 0))
        bc = np.zeros((points.shape[0], self.dimension + 1))
        tetras = np.zeros(points.shape[0], dtype="int64")
        inside = np.zeros(points.shape[0], dtype=bool)
        npts = 0
        npts_step = int(1e4)
        # break into blocks of 10k points
        while npts < points.shape[0]:
            cell_index, inside = self.aabb_grid.position_to_cell_index(
                points[: npts + npts_step, :]
            )
            global_index = self.aabb_grid.global_cell_indices(cell_index)
            tetra_indices = self.aabb_table[global_index[inside], :].tocoo()
            # tetra_indices[:] = -1
            row = tetra_indices.row
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in

            vertices = self.nodes[self.elements[col, : self.dimension + 1]]
            pos = points[row, : self.dimension]
            row = tetra_indices.row
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in
            vpa = pos[:, :] - vertices[:, 0, :]
            vba = vertices[:, 1, :] - vertices[:, 0, :]
            vca = vertices[:, 2, :] - vertices[:, 0, :]
            d00 = np.einsum('ij,ij->i', vba, vba)
            d01 = np.einsum('ij,ij->i', vba, vca)
            d11 = np.einsum('ij,ij->i', vca, vca)
            d20 = np.einsum('ij,ij->i', vpa, vba)
            d21 = np.einsum('ij,ij->i', vpa, vca)
            denom = d00 * d11 - d01 * d01
            c = np.zeros((denom.shape[0], 3))
            c[:, 0] = (d11 * d20 - d01 * d21) / denom
            c[:, 1] = (d00 * d21 - d01 * d20) / denom
            c[:, 2] = 1.0 - c[:, 0] - c[:, 1]

            mask = np.all(c >= 0, axis=1)
            if return_verts:
                verts[: npts + npts_step, :, :][row[mask], :, :] = vertices[mask, :, :]
            bc[: npts + npts_step, :][row[mask], :] = c[mask, :]
            tetras[: npts + npts_step][row[mask]] = col[mask]
            inside[: npts + npts_step][row[mask]] = True
            npts += npts_step
        tetra_return = np.zeros((points.shape[0])).astype(int)
        tetra_return[:] = -1
        tetra_return[inside] = tetras[inside]
        return verts, bc, tetra_return, inside

    def get_element_gradient_for_location(
        self, pos: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get the element gradients for a location

        Parameters
        ----------
        pos : np.array
            location to evaluate

        Returns
        -------

        """
        verts, c, tri, inside = self.get_element_for_location(pos, return_verts=False)
        return self.evaluate_shape_derivatives(pos, tri)

    def vtk(self, node_properties={}, cell_properties={}):
        """
        Create a vtk unstructured grid from the mesh
        """
        import pyvista as pv

        grid = pv.UnstructuredGrid()
        grid.points = self.nodes
        grid.cell_types = np.ones(self.elements.shape[0]) * pv.vtk.VTK_TRIANGLE
        grid.cells = np.c_[np.ones(self.elements.shape[0]) * 3, self.elements]
        for key, value in node_properties.items():
            grid.point_data[key] = value
        for key, value in cell_properties.items():
            grid.cell_data[key] = value
        return grid
