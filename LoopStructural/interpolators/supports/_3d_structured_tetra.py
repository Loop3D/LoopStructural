"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

import numpy as np
from ._3d_base_structured import BaseStructuredSupport
from . import SupportType
from scipy.sparse import coo_matrix, tril
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class TetMesh(BaseStructuredSupport):
    """ """

    def __init__(self, origin=np.zeros(3), nsteps=np.ones(3) * 10, step_vector=np.ones(3)):
        BaseStructuredSupport.__init__(self, origin, nsteps, step_vector)
        self.type = SupportType.TetMesh
        self.tetra_mask_even = np.array(
            [[7, 1, 2, 4], [6, 2, 4, 7], [5, 1, 4, 7], [0, 1, 2, 4], [3, 1, 2, 7]]
        )

        self.tetra_mask = np.array(
            [[0, 6, 5, 3], [7, 3, 5, 6], [4, 0, 5, 6], [2, 0, 3, 6], [1, 0, 3, 5]]
        )
        self.shared_element_relationships = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int
        )
        self.shared_elements = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int
        )
        self.cg = None
        self._elements = None

        self._init_face_table()

    def onGeometryChange(self):
        self._elements = None
        self.shared_element_relationships = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int
        )
        self.shared_elements = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int
        )
        self._init_face_table()
        if self.interpolator is not None:
            self.interpolator.reset()

    @property
    def neighbours(self):
        return self.get_neighbours()

    @property
    def ntetra(self) -> int:
        return np.prod(self.nsteps_cells) * 5

    @property
    def n_elements(self) -> int:
        return self.ntetra

    @property
    def n_cells(self) -> int:
        return np.prod(self.nsteps_cells)

    @property
    def elements(self):
        if self._elements is None:
            self._elements = self.get_elements()
        return self._elements

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
    
    @property
    def element_scale(self):
        size = self.element_size
        size-= np.min(size)
        size/= np.max(size)
        size+=1.
        return size

    @property
    def barycentre(self) -> np.ndarray:
        """
        Return the barycentres of all tetrahedrons or of specified tetras using
        global index

        Returns
        -------
        barycentres : numpy array
            barycentres of all tetrahedrons
        """

        tetra = self.elements
        barycentre = np.sum(self.nodes[tetra][:, :, :], axis=1) / 4.0
        return barycentre

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
        elements = self.elements
        neighbours = self.get_neighbours()
        # add array of bool to the location where there are elements for each node

        # use this to determine shared faces

        element_nodes = coo_matrix(
            (np.ones(elements.shape[0] * 4), (rows.ravel(), elements.ravel())),
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
    def shared_element_scale(self):
        return self.shared_element_size / np.mean(self.shared_element_size)

    def evaluate_value(self, pos: np.ndarray, property_array: np.ndarray) -> np.ndarray:
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
            c[inside, :] * property_array[self.elements[tetras[inside]]], axis=1
        )
        return values

    def evaluate_gradient(self, pos: np.ndarray, property_array: np.ndarray) -> np.ndarray:
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
            element_gradients[inside, :, :]
            * property_array[self.elements[tetras[inside]][:, None, :]]
        ).sum(2)
        # length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def inside(self, pos: np.ndarray):
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside

    def get_element_for_location(self, pos: np.ndarray):
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        pos = np.array(pos)
        pos = pos[:, : self.dimension]
        inside = self.inside(pos)
        # initialise array for tetrahedron vertices
        vertices = np.zeros((pos.shape[0], 5, 4, 3))
        vertices[:] = np.nan
        # get cell indexes
        cell_indexes, inside = self.position_to_cell_index(pos)
        # determine if using +ve or -ve mask
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        # get cell corners
        corner_indexes = self.cell_corner_indexes(cell_indexes)  # global_index_to_node_index(gi)
        # convert to node locations
        nodes = self.node_indexes_to_position(corner_indexes)

        vertices[even_mask, :, :, :] = nodes[even_mask, :, :][:, self.tetra_mask_even, :]
        vertices[~even_mask, :, :, :] = nodes[~even_mask, :, :][:, self.tetra_mask, :]
        # changing order to points, tetra, nodes, coord
        # vertices = vertices.swapaxes(0, 2)
        # vertices = vertices.swapaxes(1, 2)
        # use scalar triple product to calculate barycentric coords

        vap = pos[:, None, :] - vertices[:, :, 0, :]
        vbp = pos[:, None, :] - vertices[:, :, 1, :]
        #         # vcp = p - points[:, 2, :]
        #         # vdp = p - points[:, 3, :]
        vab = vertices[:, :, 1, :] - vertices[:, :, 0, :]
        vac = vertices[:, :, 2, :] - vertices[:, :, 0, :]
        vad = vertices[:, :, 3, :] - vertices[:, :, 0, :]
        vbc = vertices[:, :, 2, :] - vertices[:, :, 1, :]
        vbd = vertices[:, :, 3, :] - vertices[:, :, 1, :]
        va = np.einsum("ikj, ikj->ik", vbp, np.cross(vbd, vbc, axisa=2, axisb=2)) / 6.0
        vb = np.einsum("ikj, ikj->ik", vap, np.cross(vac, vad, axisa=2, axisb=2)) / 6.0
        vc = np.einsum("ikj, ikj->ik", vap, np.cross(vad, vab, axisa=2, axisb=2)) / 6.0
        vd = np.einsum("ikj, ikj->ik", vap, np.cross(vab, vac, axisa=2, axisb=2)) / 6.0
        v = np.einsum("ikj, ikj->ik", vab, np.cross(vac, vad, axisa=2, axisb=2)) / 6.0
        c = np.zeros((va.shape[0], va.shape[1], 4))
        c[:, :, 0] = va / v
        c[:, :, 1] = vb / v
        c[:, :, 2] = vc / v
        c[:, :, 3] = vd / v

        # if all coords are +ve then point is inside cell
        mask = np.all(c >= 0, axis=2)
        i, j = np.where(mask)
        ## find any cases where the point belongs to two cells
        ## just use the second cell
        pairs = dict(zip(i, j))
        mask[:] = False
        mask[list(pairs.keys()), list(pairs.values())] = True

        inside = np.logical_and(inside, np.any(mask, axis=1))
        # get cell corners
        # create mask to see which cells are even
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        # create global node index list
        gi = self.global_node_indices(corner_indexes)
        # gi = xi + yi * self.nsteps[0] + zi * self.nsteps[0] * self.nsteps[1]
        # container for tetras
        tetras = np.zeros((corner_indexes.shape[0], 5, 4)).astype(int)
        tetras[even_mask, :, :] = gi[even_mask, :][:, self.tetra_mask_even]
        tetras[~even_mask, :, :] = gi[~even_mask, :][:, self.tetra_mask]
        inside = np.logical_and(inside, self.inside(pos))
        vertices_return = np.zeros((pos.shape[0], 4, 3))
        vertices_return[:] = np.nan
        # set all masks not inside to False
        mask[~inside, :] = False
        vertices_return[inside, :, :] = vertices[mask, :, :]  # [mask,:,:]#[inside,:,:]
        c_return = np.zeros((pos.shape[0], 4))
        c_return[:] = np.nan
        c_return[inside] = c[mask]
        tetra_return = np.zeros((pos.shape[0])).astype(int)
        tetra_return[:] = -1
        local_tetra_index = np.tile(np.arange(0, 5)[None, :], (mask.shape[0], 1))
        local_tetra_index = local_tetra_index[mask]
        tetra_global_index = self.tetra_global_index(cell_indexes[inside, :], local_tetra_index)
        tetra_return[inside] = tetra_global_index
        return vertices_return, c_return, tetra_return, inside

    def evaluate_shape(self, locations):
        """
        Convenience function returning barycentric coords

        """
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        return c, elements, inside

    def get_elements(self):
        """
        Get a numpy array of all of the elements in the mesh

        Returns
        -------
        numpy array elements

        """

        x = np.arange(0, self.nsteps_cells[0])
        y = np.arange(0, self.nsteps_cells[1])
        z = np.arange(0, self.nsteps_cells[2])
        ## reverse x and z so that indexing is
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        cell_indexes = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        # get cell corners
        cell_corners = self.cell_corner_indexes(cell_indexes)
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        gi = self.global_node_indices(cell_corners)
        tetras = np.zeros((cell_corners.shape[0], 5, 4)).astype("int64")
        tetras[even_mask, :, :] = gi[even_mask, :][:, self.tetra_mask_even]
        tetras[~even_mask, :, :] = gi[~even_mask, :][:, self.tetra_mask]

        return tetras.reshape((tetras.shape[0] * tetras.shape[1], tetras.shape[2]))

    def tetra_global_index(self, indices, tetra_index):
        """
        Get the global index of a tetra from the cell index and the local tetra index

        Parameters
        ----------
        indices
        tetra_index

        Returns
        -------

        """
        return (
            tetra_index
            + indices[:, 0] * 5
            + self.nsteps_cells[0] * indices[:, 1] * 5
            + self.nsteps_cells[0] * self.nsteps_cells[1] * indices[:, 2] * 5
        )

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
            elements = np.arange(0, self.ntetra)
        x = np.arange(0, self.nsteps_cells[0])
        y = np.arange(0, self.nsteps_cells[1])
        z = np.arange(0, self.nsteps_cells[2])

        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        cell_indexes = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        # c_xi = c_xi.flatten(order="F")
        # c_yi = c_yi.flatten(order="F")
        # c_zi = c_zi.flatten(order="F")
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        # get cell corners
        corner_indexes = self.cell_corner_indexes(cell_indexes)  # global_index_to_node_index(gi)
        # convert to node locations
        nodes = self.node_indexes_to_position(corner_indexes)

        points = np.zeros((self.n_cells, 5, 4, 3))
        points[even_mask, :, :, :] = nodes[even_mask, :, :][:, self.tetra_mask_even, :]
        points[~even_mask, :, :, :] = nodes[~even_mask, :, :][:, self.tetra_mask, :]

        # changing order to points, tetra, nodes, coord
        # points = points.swapaxes(0, 2)
        # points = points.swapaxes(1, 2)

        ps = points.reshape(points.shape[0] * points.shape[1], points.shape[2], points.shape[3])

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

    def evaluate_shape_derivatives(self, pos, elements=None):
        inside = None
        if elements is not None:
            inside = np.ones(elements.shape[0], dtype=bool)
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(pos)
            # np.arange(0, self.n_elements, dtype=int)

        return (
            self.get_element_gradients(elements),
            elements,
            inside,
        )

    def get_element_gradient_for_location(self, pos: np.ndarray):
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
        # m[~inside,:,:] = np.nan
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.zeros_like(m)
        element_gradients[:] = np.nan
        element_gradients[inside, :, :] = np.linalg.inv(m[inside, :, :])
        # element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        return vertices, element_gradients, tetras, inside

    def global_node_indicies(self, indexes: np.ndarray):
        """
        Convert from node indexes to global node index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        indexes = np.array(indexes).swapaxes(0, 2)
        return (
            indexes[:, :, 0]
            + self.nsteps[None, None, 0] * indexes[:, :, 1]
            + self.nsteps[None, None, 0] * self.nsteps[None, None, 1] * indexes[:, :, 2]
        )

    def global_cell_indicies(self, indexes: np.ndarray):
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        indexes = np.array(indexes).swapaxes(0, 2)
        return (
            indexes[:, :, 0]
            + self.nsteps_cells[None, None, 0] * indexes[:, :, 1]
            + self.nsteps_cells[None, None, 0] * self.nsteps_cells[None, None, 1] * indexes[:, :, 2]
        )

    def global_index_to_node_index(self, global_index: np.ndarray):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        # determine the ijk indices for the global index.
        # remainder when dividing by nx = i
        # remained when dividing modulus of nx by ny is j
        x_index = global_index % self.nsteps[0, None]
        y_index = global_index // self.nsteps[0, None] % self.nsteps[1, None]
        z_index = global_index // self.nsteps[0, None] // self.nsteps[1, None]
        return x_index, y_index, z_index

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        # determine the ijk indices for the global index.
        # remainder when dividing by nx = i
        # remained when dividing modulus of nx by ny is j

        x_index = global_index % self.nsteps_cells[0, None]
        y_index = global_index // self.nsteps_cells[0, None] % self.nsteps_cells[1, None]
        z_index = global_index // self.nsteps_cells[0, None] // self.nsteps_cells[1, None]
        return x_index, y_index, z_index

    def get_neighbours(self) -> np.ndarray:
        """
        This function goes through all of the elements in the mesh and assembles a numpy array
        with the neighbours for each element

        Returns
        -------

        """
        # elements = self.get_elements()
        # neighbours = np.zeros((self.ntetra,4)).astype('int64')
        # neighbours[:] = -1
        # tetra_neighbours(elements,neighbours)
        # return neighbours
        tetra_index = np.arange(0, self.ntetra)
        neighbours = np.zeros((self.ntetra, 4)).astype("int64")
        neighbours[:] = -9999
        neighbours[tetra_index % 5 == 0, :] = (
            tetra_index[tetra_index % 5 == 0, None] + np.arange(1, 5)[None, :]
        )  # first tetra is the centre one so all of its neighbours are in the same cell
        neighbours[tetra_index % 5 != 0, 0] = np.tile(
            tetra_index[tetra_index % 5 == 0], (4, 1)
        ).flatten(
            order="F"
        )  # add first tetra to other neighbours

        # now create masks for the different tetra indexes
        one_mask = tetra_index % 5 == 1
        two_mask = tetra_index % 5 == 2
        three_mask = tetra_index % 5 == 3
        four_mask = tetra_index % 5 == 4

        # create masks for whether cell is odd or even
        odd_mask = np.sum(self.global_index_to_cell_index(tetra_index // 5), axis=0) % 2 == 1
        odd_mask = odd_mask.astype(bool)

        # apply masks to
        masks = []
        masks.append(
            [
                np.logical_and(one_mask, odd_mask),
                np.array([[1, 0, 0, 1], [0, 1, 0, 2], [0, 0, 1, 4]]),
            ]
        )
        masks.append(
            [
                np.logical_and(two_mask, odd_mask),
                np.array([[-1, 0, 0, 2], [0, -1, 0, 1], [0, 0, 1, 3]]),
            ]
        )
        masks.append(
            [
                np.logical_and(three_mask, odd_mask),
                np.array([[-1, 0, 0, 4], [0, 1, 0, 3], [0, 0, -1, 1]]),
            ]
        )
        masks.append(
            [
                np.logical_and(four_mask, odd_mask),
                np.array([[1, 0, 0, 3], [0, -1, 0, 4], [0, 0, -1, 2]]),
            ]
        )

        masks.append(
            [
                np.logical_and(one_mask, ~odd_mask),
                np.array([[-1, 0, 0, 1], [0, 1, 0, 2], [0, 0, 1, 3]]),
            ]
        )
        masks.append(
            [
                np.logical_and(two_mask, ~odd_mask),
                np.array([[1, 0, 0, 2], [0, -1, 0, 1], [0, 0, 1, 4]]),
            ]
        )
        masks.append(
            [
                np.logical_and(three_mask, ~odd_mask),
                np.array([[-1, 0, 0, 4], [0, -1, 0, 3], [0, 0, -1, 2]]),
            ]
        )
        masks.append(
            [
                np.logical_and(four_mask, ~odd_mask),
                np.array([[1, 0, 0, 3], [0, 1, 0, 4], [0, 0, -1, 1]]),
            ]
        )

        for m in masks:
            logic = m[0]
            mask = m[1]
            c_xi, c_yi, c_zi = self.global_index_to_cell_index(tetra_index[logic] // 5)
            # mask = np.array([[1,0,0,4],[0,0,-1,2],[0,1,0,3],[0,0,0,0]])
            neigh_cell = np.zeros((c_xi.shape[0], 3, 3)).astype(int)
            neigh_cell[:, :, 0] = c_xi[:, None] + mask[:, 0]
            neigh_cell[:, :, 1] = c_yi[:, None] + mask[:, 1]
            neigh_cell[:, :, 2] = c_zi[:, None] + mask[:, 2]
            inside = neigh_cell[:, :, 0] >= 0
            inside = np.logical_and(inside, neigh_cell[:, :, 1] >= 0)
            inside = np.logical_and(inside, neigh_cell[:, :, 2] >= 0)
            inside = np.logical_and(inside, neigh_cell[:, :, 0] < self.nsteps_cells[0])
            inside = np.logical_and(inside, neigh_cell[:, :, 1] < self.nsteps_cells[1])
            inside = np.logical_and(inside, neigh_cell[:, :, 2] < self.nsteps_cells[2])

            global_neighbour_idx = np.zeros((c_xi.shape[0], 4)).astype(int)
            global_neighbour_idx[:] = -1
            global_neighbour_idx = (
                neigh_cell[:, :, 0]
                + neigh_cell[:, :, 1] * self.nsteps_cells[0]
                + neigh_cell[:, :, 2] * self.nsteps_cells[0] * self.nsteps_cells[1]
            ) * 5 + mask[:, 3]
            global_neighbour_idx[~inside] = -1
            neighbours[logic, 1:] = global_neighbour_idx

        return neighbours

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
        for prop in node_properties:
            grid[prop] = node_properties[prop]
        for prop in cell_properties:
            grid.cell_arrays[prop] = cell_properties[prop]
        return grid
