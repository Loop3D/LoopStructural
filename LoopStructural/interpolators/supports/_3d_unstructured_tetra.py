"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

from ast import Tuple


import numpy as np
from scipy.sparse import csr_matrix, coo_matrix, tril

from . import StructuredGrid
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
        self._nodes = np.array(nodes)
        if self._nodes.shape[1] != 3:
            raise ValueError("Nodes must be 3D")
        self.neighbours = np.array(neighbours, dtype=np.int64)
        if self.neighbours.shape[1] != 4:
            raise ValueError("Neighbours array is too big")
        self._elements = np.array(elements, dtype=np.int64)
        if self.elements.shape[0] != self.neighbours.shape[0]:
            raise ValueError("Number of elements and neighbours do not match")
        self._barycentre = np.sum(self.nodes[self.elements[:, :4]][:, :, :], axis=1) / 4.0
        self.minimum = np.min(self.nodes, axis=0)
        self.maximum = np.max(self.nodes, axis=0)
        length = self.maximum - self.minimum
        self.minimum -= length * 0.1
        self.maximum += length * 0.1
        if self.elements.shape[0] < 2000:
            self.aabb_grid = StructuredGrid(self.minimum, nsteps=[2, 2, 2], step_vector=[1, 1, 1])
        else:
            if aabb_nsteps is None:
                box_vol = np.prod(self.maximum - self.minimum)
                element_volume = box_vol / (len(self.elements) / 20)
                # calculate the step vector of a regular cube
                step_vector = np.zeros(3)
                step_vector[:] = element_volume ** (1.0 / 3.0)
                # number of steps is the length of the box / step vector
                aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(int)
                # make sure there is at least one cell in every dimension
                aabb_nsteps[aabb_nsteps < 2] = 2
            aabb_nsteps = np.array(aabb_nsteps, dtype=int)
            step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
            self.aabb_grid = StructuredGrid(
                self.minimum, nsteps=aabb_nsteps, step_vector=step_vector
            )
        # make a big table to store which tetra are in which element.
        # if this takes up too much memory it could be simplified by using sparse matrices or dict but
        # at the expense of speed
        self.aabb_table = csr_matrix((self.aabb_grid.n_elements, len(self.elements)), dtype=bool)
        self.shared_element_relationships = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int
        )
        self.shared_elements = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int
        )
        self._init_face_table()
        self._initialise_aabb()

    def set_nelements(self, nelements):
        raise NotImplementedError("Cannot set number of elements for unstructured mesh")

    @property
    def nodes(self):
        return self._nodes

    @property
    def elements(self):
        return self._elements

    @property
    def barycentre(self):
        return self._barycentre

    @property
    def n_nodes(self):
        return self.nodes.shape[0]

    def onGeometryChange(self):
        pass

    def _init_face_table(self):
        """
        Fill table containing elements that share a face, and another
        table that contains the nodes for a face.
        """
        # need to identify the shared nodes for pairs of elements
        # we do this by creating a sparse matrix that has N rows (number of elements)
        # and M columns (number of nodes).
        # We then fill the location where a node is in an element with true
        # Then we create a table for the pairs of elements in the mesh
        # we have the neighbour relationships, which are the 4 neighbours for each element
        # create a new table that shows the element index repeated four times
        # flatten both of these arrays so we effectively have a table with pairs of neighbours
        # disgard the negative neighbours because these are border neighbours
        rows = np.tile(np.arange(self.n_elements)[:, None], (1, 4))
        elements = self.get_elements()
        neighbours = self.get_neighbours()
        # add array of bool to the location where there are elements for each node

        # use this to determine shared faces

        element_nodes = coo_matrix(
            (np.ones(elements.shape[0] * 4), (rows.ravel(), elements[:, :4].ravel())),
            shape=(self.n_elements, self.n_nodes),
            dtype=bool,
        ).tocsr()
        n1 = np.tile(np.arange(neighbours.shape[0], dtype=int)[:, None], (1, 4))
        n1 = n1.flatten()
        n2 = neighbours.flatten()
        n1 = n1[n2 >= 0]
        n2 = n2[n2 >= 0]
        el_rel = np.zeros((self.neighbours.flatten().shape[0], 2), dtype=int)
        el_rel[:] = -1
        el_rel[np.arange(n1.shape[0]), 0] = n1
        el_rel[np.arange(n1.shape[0]), 1] = n2
        el_rel = el_rel[el_rel[:, 0] >= 0, :]

        # el_rel2 = np.zeros((self.neighbours.flatten().shape[0], 2), dtype=int)
        self.shared_element_relationships[:] = -1
        el_pairs = coo_matrix((np.ones(el_rel.shape[0]), (el_rel[:, 0], el_rel[:, 1]))).tocsr()
        i, j = tril(el_pairs).nonzero()
        self.shared_element_relationships[: len(i), 0] = i
        self.shared_element_relationships[: len(i), 1] = j

        self.shared_element_relationships = self.shared_element_relationships[
            self.shared_element_relationships[:, 0] >= 0, :
        ]

        faces = element_nodes[self.shared_element_relationships[:, 0], :].multiply(
            element_nodes[self.shared_element_relationships[:, 1], :]
        )
        shared_faces = faces[np.array(np.sum(faces, axis=1) == 3).flatten(), :]
        row, col = shared_faces.nonzero()
        row = row[row.argsort()]
        col = col[row.argsort()]
        shared_face_index = np.zeros((shared_faces.shape[0], 3), dtype=int)
        shared_face_index[:] = -1
        shared_face_index[row.reshape(-1, 3)[:, 0], :] = col.reshape(-1, 3)

        self.shared_elements[np.arange(self.shared_element_relationships.shape[0]), :] = (
            shared_face_index
        )
        # resize
        self.shared_elements = self.shared_elements[: len(self.shared_element_relationships), :]
        # flag = np.zeros(self.elements.shape[0])
        # face_index = 0
        # for i, t in enumerate(self.elements):
        #     flag[i] = True
        #     for n in self.neighbours[i]:
        #         if n < 0:
        #             continue
        #         if flag[n]:
        #             continue
        #         face_node_index = 0
        #         self.shared_element_relationships[face_index, 0] = i
        #         self.shared_element_relationships[face_index, 1] = n
        #         for v in t:
        #             if v in self.elements[n, :4]:
        #                 self.shared_elements[face_index, face_node_index] = v
        #                 face_node_index += 1

        #         face_index += 1
        # self.shared_elements = self.shared_elements[:face_index, :]
        # self.shared_element_relationships = self.shared_element_relationships[
        #     :face_index, :
        # ]

    def _initialise_aabb(self):
        """assigns the tetras to the grid cells where the bounding box
        of the tetra element overlaps the grid cell.
        It could be changed to use the separating axis theorem, however this would require
        significantly more calculations. (12 more I think).. #TODO test timing
        """
        # calculate the bounding box for all tetraherdon in the mesh
        # find the min/max extents for xyz
        # tetra_bb = np.zeros((self.elements.shape[0], 19, 3))
        minx = np.min(self.nodes[self.elements[:, :4], 0], axis=1)
        maxx = np.max(self.nodes[self.elements[:, :4], 0], axis=1)
        miny = np.min(self.nodes[self.elements[:, :4], 1], axis=1)
        maxy = np.max(self.nodes[self.elements[:, :4], 1], axis=1)
        minz = np.min(self.nodes[self.elements[:, :4], 2], axis=1)
        maxz = np.max(self.nodes[self.elements[:, :4], 2], axis=1)
        cell_indexes = self.aabb_grid.global_index_to_cell_index(
            np.arange(self.aabb_grid.n_elements)
        )
        corners = self.aabb_grid.cell_corner_indexes(cell_indexes)
        positions = self.aabb_grid.node_indexes_to_position(corners)
        ## Because we known the node orders just select min/max from each
        # coordinate. Use these to check whether the tetra is in the cell
        x_boundary = positions[:, [0, 1], 0]
        y_boundary = positions[:, [0, 2], 1]
        z_boundary = positions[:, [0, 6], 2]
        a = np.logical_and(
            minx[None, :] > x_boundary[:, None, 0],
            minx[None, :] < x_boundary[:, None, 1],
        )  # min point between cell
        b = np.logical_and(
            maxx[None, :] < x_boundary[:, None, 1],
            maxx[None, :] > x_boundary[:, None, 0],
        )  # max point between cell
        c = np.logical_and(
            minx[None, :] < x_boundary[:, None, 0],
            maxx[None, :] > x_boundary[:, None, 0],
        )  # min point < than cell & max point > cell

        x_logic = np.logical_or(np.logical_or(a, b), c)

        a = np.logical_and(
            miny[None, :] > y_boundary[:, None, 0],
            miny[None, :] < y_boundary[:, None, 1],
        )  # min point between cell
        b = np.logical_and(
            maxy[None, :] < y_boundary[:, None, 1],
            maxy[None, :] > y_boundary[:, None, 0],
        )  # max point between cell
        c = np.logical_and(
            miny[None, :] < y_boundary[:, None, 0],
            maxy[None, :] > y_boundary[:, None, 0],
        )  # min point < than cell & max point > cell

        y_logic = np.logical_or(np.logical_or(a, b), c)

        a = np.logical_and(
            minz[None, :] > z_boundary[:, None, 0],
            minz[None, :] < z_boundary[:, None, 1],
        )  # min point between cell
        b = np.logical_and(
            maxz[None, :] < z_boundary[:, None, 1],
            maxz[None, :] > z_boundary[:, None, 0],
        )  # max point between cell
        c = np.logical_and(
            minz[None, :] < z_boundary[:, None, 0],
            maxz[None, :] > z_boundary[:, None, 0],
        )  # min point < than cell & max point > cell

        z_logic = np.logical_or(np.logical_or(a, b), c)
        logic = np.logical_and(x_logic, y_logic)
        logic = np.logical_and(logic, z_logic)

        self.aabb_table = csr_matrix(logic)

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
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        v2 = self.nodes[elements[:, 2], :] - self.nodes[elements[:, 0], :]
        return np.cross(v1, v2, axisa=1, axisb=1)

    @property
    def shared_element_size(self):
        """
        Get the area of the share triangle
        """
        norm = self.shared_element_norm
        return 0.5 * np.linalg.norm(norm, axis=1)

    @property
    def element_size(self):
        """Calculate the volume of a tetrahedron using the 4 corners
        volume = abs(det(A))/6 where A is the jacobian of the corners

        Returns
        -------
        _type_
            _description_
        """
        vecs = (
            self.nodes[self.elements[:, :4], :][:, 1:, :]
            - self.nodes[self.elements[:, :4], :][:, 0, None, :]
        )
        return np.abs(np.linalg.det(vecs)) / 6

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

        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside

    def get_elements(self):
        return self.elements

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
            cell_index, inside = self.aabb_grid.position_to_cell_index(
                points[: npts + npts_step, :]
            )
            global_index = (
                cell_index[:, 0]
                + self.aabb_grid.nsteps_cells[None, 0] * cell_index[:, 1]
                + self.aabb_grid.nsteps_cells[None, 0]
                * self.aabb_grid.nsteps_cells[None, 1]
                * cell_index[:, 2]
            )

            tetra_indices = self.aabb_table[global_index[inside], :].tocoo()
            # tetra_indices[:] = -1
            row = tetra_indices.row
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in
            vertices = self.nodes[self.elements[col, :4]]
            pos = points[row, :]
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

            verts[: npts + npts_step, :, :][row[mask], :, :] = vertices[mask, :, :]
            bc[: npts + npts_step, :][row[mask], :] = c[mask, :]
            tetras[: npts + npts_step][row[mask]] = col[mask]
            inside[: npts + npts_step][row[mask]] = True
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
        # points = np.zeros((5, 4, self.n_cells, 3))
        # points[:, :, even_mask, :] = nodes[:, even_mask, :][self.tetra_mask_even, :, :]
        # points[:, :, ~even_mask, :] = nodes[:, ~even_mask, :][self.tetra_mask, :, :]

        # # changing order to points, tetra, nodes, coord
        # points = points.swapaxes(0, 2)
        # points = points.swapaxes(1, 2)
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
        return self.neighbours

    def vtk(self, node_properties={}, cell_properties={}):
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
