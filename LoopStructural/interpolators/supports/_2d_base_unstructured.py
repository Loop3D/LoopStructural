"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

from abc import abstractmethod
import logging
from typing import Tuple
import numpy as np

from LoopStructural.geometry import UnstructuredMesh2DGeometry
from . import SupportType
from ._base_support import BaseSupport

logger = logging.getLogger(__name__)


class BaseUnstructured2d(BaseSupport):
    """ """

    dimension = 2

    def __init__(self, elements, vertices, neighbours, aabb_nsteps=None):
        self.type = SupportType.BaseUnstructured2d
        self._geom = UnstructuredMesh2DGeometry(
            elements, vertices, neighbours, aabb_nsteps=aabb_nsteps
        )
        if self.elements.shape[1] == 3:
            self.order = 1
        elif self.elements.shape[1] == 6:
            self.order = 2
        self.dof = self.vertices.shape[0]

    @property
    def vertices(self):
        return self._geom.vertices

    @property
    def neighbours(self):
        return self._geom.neighbours

    @property
    def minimum(self):
        return self._geom.minimum

    @property
    def maximum(self):
        return self._geom.maximum

    @property
    def aabb_grid(self):
        return self._geom.aabb_grid

    @property
    def aabb_table(self):
        return self._geom.aabb_table

    def set_nelements(self, nelements) -> int:
        raise NotImplementedError

    @property
    def shared_elements(self):
        return self._geom.shared_elements

    @property
    def shared_element_relationships(self):
        return self._geom.shared_element_relationships

    @property
    def elements(self):
        return self._geom.elements

    def onGeometryChange(self):
        pass

    @property
    def n_elements(self):
        return self._geom.n_elements

    @property
    def n_nodes(self):
        return self._geom.n_nodes

    def inside(self, pos):
        if pos.shape[1] > self.dimension:
            logger.warning(f"Converting {pos.shape[1]} to 3d using first {self.dimension} columns")
            pos = pos[:, : self.dimension]
        return self._geom.inside(pos)

    @property
    def ncps(self):
        """
        Returns the number of nodes for an element in the mesh
        """
        return self._geom.ncps

    @property
    def nodes(self):
        """
        Gets the nodes of the mesh as a property rather than using a function, accessible as a property! Python magic!

        Returns
        -------
        nodes : np.array((N,3))
            Fortran ordered
        """
        return self._geom.nodes

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
        return self._geom.barycentre

    @property
    def shared_element_norm(self):
        """
        Get the normal to all of the shared elements
        """
        return self._geom.shared_element_norm

    @property
    def shared_element_size(self):
        """
        Get the size of the shared elements
        """
        return self._geom.shared_element_size

    @property
    def element_size(self):
        return self._geom.element_size

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
        return self._geom.element_area(elements)

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
            chunk = points[npts : npts + npts_step, :]
            cell_index, chunk_inside = self.aabb_grid.position_to_cell_index(chunk)
            global_index = self.aabb_grid.global_cell_indices(cell_index)
            tetra_indices = self.aabb_table[global_index[chunk_inside], :].tocoo()
            # tetra_indices[:] = -1
            row = tetra_indices.row
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in

            vertices = self.nodes[self.elements[col, : self.dimension + 1]]
            pos = chunk[row, : self.dimension]
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
            # d11*d20-d01*d21 and d00*d21-d01*d20 are the barycentric weights
            # for vertices[:,1] and vertices[:,2] respectively (standard
            # Ericson barycentric technique) - assign to the matching columns
            # so that c[:, i] lines up with self.elements[tri, i] everywhere
            # else in the codebase (evaluate_value, add_value_constraints, ...)
            c[:, 1] = (d11 * d20 - d01 * d21) / denom
            c[:, 2] = (d00 * d21 - d01 * d20) / denom
            c[:, 0] = 1.0 - c[:, 1] - c[:, 2]

            mask = np.all(c >= 0, axis=1)
            if return_verts:
                verts[npts : npts + npts_step, :, :][row[mask], :, :] = vertices[mask, :, :]
            bc[npts : npts + npts_step, :][row[mask], :] = c[mask, :]
            tetras[npts : npts + npts_step][row[mask]] = col[mask]
            inside[npts : npts + npts_step][row[mask]] = True
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

    def vtk(self, node_properties=None, cell_properties=None):
        """
        Create a vtk unstructured grid from the mesh
        """
        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
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
