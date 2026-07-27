"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

from ast import Tuple


import numpy as np

from LoopStructural.geometry import UnstructuredMeshGeometry
from LoopStructural.utils import getLogger
from . import SupportType
from ._base_support import BaseSupport

logger = getLogger(__name__)


class UnStructuredTetMesh(BaseSupport):
    """ """

    dimension = 3

    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        neighbours: np.ndarray,
        aabb_nsteps=None,
    ):
        """An unstructured mesh defined by nodes, elements and neighbours
        An axis aligned bounding box (AABB) is used to speed up finding
        which tetra a point is in.
        The aabb grid is calculated so that there are approximately 10 tetra per
        element.

        Parameters
        ----------
        nodes : array or array like
            container of vertex locations
        elements : array or array like, dtype cast to long
            container of tetra indicies
        neighbours : array or array like, dtype cast to long
            array containing element neighbours
        aabb_nsteps : list, optional
            force nsteps for aabb, by default None
        """
        self.type = SupportType.UnStructuredTetMesh
        self._geom = UnstructuredMeshGeometry(nodes, elements, neighbours, aabb_nsteps=aabb_nsteps)

    def set_nelements(self, nelements):
        raise NotImplementedError("Cannot set number of elements for unstructured mesh")

    @property
    def nodes(self):
        return self._geom.nodes

    @property
    def elements(self):
        return self._geom.elements

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

    @property
    def shared_elements(self):
        return self._geom.shared_elements

    @property
    def shared_element_relationships(self):
        return self._geom.shared_element_relationships

    @property
    def barycentre(self):
        return self._geom.barycentre

    @property
    def n_nodes(self):
        return self._geom.n_nodes

    def onGeometryChange(self):
        pass

    @property
    def ntetra(self):
        return self.elements.shape[0]

    @property
    def n_elements(self):
        return self.ntetra

    @property
    def n_cells(self):
        return None

    @property
    def shared_element_norm(self):
        """
        Get the normal to all of the shared elements
        """
        return self._geom.shared_element_norm

    @property
    def shared_element_size(self):
        """
        Get the area of the share triangle
        """
        return self._geom.shared_element_size

    @property
    def element_size(self):
        """Calculate the volume of a tetrahedron using the 4 corners
        volume = abs(det(A))/6 where A is the jacobian of the corners

        Returns
        -------
        np.ndarray
            array of length n_elements containing the volume of each tetrahedron
        """
        return self._geom.element_size

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        Get the gradients of all tetras

        Parameters
        ----------
        elements

        Returns
        -------

        """
        inside = None
        if elements is not None:
            inside = np.zeros(self.n_elements, dtype=bool)
            inside[elements] = True
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(locations)
            # elements = np.arange(0, self.n_elements, dtype=int)
        ps = self.nodes[self.elements, :]
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I

        return element_gradients[elements, :, :], elements, inside

    def evaluate_shape(self, locations):
        """
        Convenience function returning barycentric coords

        """
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        return c, elements, inside

    def evaluate_value(self, pos, property_array):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - string
            property name

        Returns
        -------

        """
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        vertices, c, tetras, inside = self.get_element_for_location(pos)
        values[inside] = np.sum(
            c[inside, :] * property_array[self.elements[tetras[inside], :]], axis=1
        )
        return values

    def evaluate_gradient(self, pos, property_array):
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
        values = np.zeros(pos.shape)
        values[:] = np.nan
        (
            vertices,
            element_gradients,
            tetras,
            inside,
        ) = self.get_element_gradient_for_location(pos)
        # grads = np.zeros(tetras.shape)
        values[inside, :] = (
            element_gradients[inside, :, :] * property_array[self.elements[tetras][inside, None, :]]
        ).sum(2)
        # length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def inside(self, pos):
        if pos.shape[1] > 3:
            logger.warning(f"Converting {pos.shape[1]} to 3d using first 3 columns")
            pos = pos[:, :3]
        return self._geom.inside(pos)

    def get_elements(self):
        return self._geom.get_elements()

    def get_element_for_location(self, points: np.ndarray) -> Tuple:
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        verts = np.zeros((points.shape[0], 4, 3))
        bc = np.zeros((points.shape[0], 4))
        tetras = np.zeros(points.shape[0], dtype="int64")
        inside = np.zeros(points.shape[0], dtype=bool)
        npts = 0
        npts_step = int(1e4)
        # break into blocks of 10k points
        while npts < points.shape[0]:
            chunk = points[npts : npts + npts_step, :]
            cell_index, chunk_inside = self.aabb_grid.position_to_cell_index(chunk)
            global_index = (
                cell_index[:, 0]
                + self.aabb_grid.nsteps_cells[None, 0] * cell_index[:, 1]
                + self.aabb_grid.nsteps_cells[None, 0]
                * self.aabb_grid.nsteps_cells[None, 1]
                * cell_index[:, 2]
            )

            tetra_indices = self.aabb_table[global_index[chunk_inside], :].tocoo()
            # tetra_indices[:] = -1
            row = tetra_indices.row
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in
            vertices = self.nodes[self.elements[col, :4]]
            pos = chunk[row, :]
            vap = pos[:, :] - vertices[:, 0, :]
            vbp = pos[:, :] - vertices[:, 1, :]
            #         # vcp = p - points[:, 2, :]
            #         # vdp = p - points[:, 3, :]
            vab = vertices[:, 1, :] - vertices[:, 0, :]
            vac = vertices[:, 2, :] - vertices[:, 0, :]
            vad = vertices[:, 3, :] - vertices[:, 0, :]
            vbc = vertices[:, 2, :] - vertices[:, 1, :]
            vbd = vertices[:, 3, :] - vertices[:, 1, :]

            va = np.einsum("ij, ij->i", vbp, np.cross(vbd, vbc, axisa=1, axisb=1)) / 6.0
            vb = np.einsum("ij, ij->i", vap, np.cross(vac, vad, axisa=1, axisb=1)) / 6.0
            vc = np.einsum("ij, ij->i", vap, np.cross(vad, vab, axisa=1, axisb=1)) / 6.0
            vd = np.einsum("ij, ij->i", vap, np.cross(vab, vac, axisa=1, axisb=1)) / 6.0
            v = np.einsum("ij, ij->i", vab, np.cross(vac, vad, axisa=1, axisb=1)) / 6.0
            c = np.zeros((va.shape[0], 4))
            c[:, 0] = va / v
            c[:, 1] = vb / v
            c[:, 2] = vc / v
            c[:, 3] = vd / v
            # inside = np.ones(c.shape[0],dtype=bool)
            mask = np.all(c >= 0, axis=1)

            verts[npts : npts + npts_step, :, :][row[mask], :, :] = vertices[mask, :, :]
            bc[npts : npts + npts_step, :][row[mask], :] = c[mask, :]
            tetras[npts : npts + npts_step][row[mask]] = col[mask]
            inside[npts : npts + npts_step][row[mask]] = True
            npts += npts_step
        tetra_return = np.zeros((points.shape[0])).astype(int)
        tetra_return[:] = -1

        tetra_return[inside] = tetras[inside]
        return verts, bc, tetra_return, inside

    def get_element_gradients(self, elements=None):
        """
        Get the gradients of all tetras

        Parameters
        ----------
        elements

        Returns
        -------

        """
        if elements is None:
            elements = np.arange(0, self.n_elements, dtype=int)
        ps = self.nodes[
            self.elements, :
        ]  # points.reshape(points.shape[0] * points.shape[1], points.shape[2], points.shape[3])
        # vertices = self.nodes[self.elements[col,:]]
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I

        return element_gradients[elements, :, :]

    def get_element_gradient_for_location(self, pos):
        """
        Get the gradient of the tetra for a location

        Parameters
        ----------
        pos

        Returns
        -------

        """
        vertices, bc, tetras, inside = self.get_element_for_location(pos)
        ps = vertices
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        return vertices, element_gradients, tetras, inside

    def get_neighbours(self):
        """
        This function goes through all of the elements in the mesh and assembles a numpy array
        with the neighbours for each element

        Returns
        -------

        """
        return self._geom.get_neighbours()

    def vtk(self, node_properties=None, cell_properties=None):
        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")

        from pyvista import CellType

        celltype = np.full(self.elements.shape[0], CellType.TETRA, dtype=np.uint8)
        elements = np.hstack(
            [np.zeros(self.elements.shape[0], dtype=int)[:, None] + 4, self.elements]
        )
        elements = elements.flatten()
        grid = pv.UnstructuredGrid(elements, celltype, self.nodes)
        for key, value in node_properties.items():
            grid[key] = value
        for key, value in cell_properties.items():
            grid.cell_arrays[key] = value

        return grid
