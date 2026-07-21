"""
Pure unstructured mesh geometry: nodes/elements/neighbours containers for
tetrahedral (3D) and triangular (2D) meshes, plus an axis-aligned bounding-box
(AABB) grid used to accelerate point-in-element lookups.
"""

import numpy as np
from scipy import sparse

from ._aabb import _initialise_aabb
from ._face_table import _init_face_table
from ._structured_grid_3d import StructuredGrid3DGeometry
from ._structured_grid_2d import StructuredGrid2DGeometry


class UnstructuredMeshGeometry:
    """An unstructured tetrahedral mesh defined by nodes, elements and neighbours.

    An axis aligned bounding box (AABB) is used to speed up finding
    which tetra a point is in. The aabb grid is calculated so that there
    are approximately 10 tetra per element.
    """

    dimension = 3

    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        neighbours: np.ndarray,
        aabb_nsteps=None,
    ):
        """

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
            self.aabb_grid = StructuredGrid3DGeometry(
                self.minimum, nsteps=[2, 2, 2], step_vector=[1, 1, 1]
            )
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
            self.aabb_grid = StructuredGrid3DGeometry(
                self.minimum, nsteps=aabb_nsteps, step_vector=step_vector
            )
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
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int
        )

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

    @property
    def n_elements(self):
        return self.elements.shape[0]

    @property
    def aabb_table(self):
        if np.sum(self._aabb_table) == 0:
            _initialise_aabb(self)
        return self._aabb_table

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

    def get_elements(self):
        return self.elements

    def get_neighbours(self):
        """
        This function goes through all of the elements in the mesh and assembles a numpy array
        with the neighbours for each element

        Returns
        -------

        """
        return self.neighbours

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
        np.ndarray
            array of length n_elements containing the volume of each tetrahedron
        """
        vecs = (
            self.nodes[self.elements[:, :4], :][:, 1:, :]
            - self.nodes[self.elements[:, :4], :][:, 0, None, :]
        )
        return np.abs(np.linalg.det(vecs)) / 6

    def inside(self, pos):
        if pos.shape[1] > 3:
            pos = pos[:, :3]

        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.minimum[None, i]
            inside *= pos[:, i] < self.maximum[None, i]
        return inside


class UnstructuredMesh2DGeometry:
    """An unstructured triangular mesh defined by vertices, elements and neighbours.

    An axis aligned bounding box (AABB) is used to speed up finding
    which triangle a point is in.
    """

    dimension = 2

    def __init__(self, elements, vertices, neighbours, aabb_nsteps=None):
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
        self.aabb_grid = StructuredGrid2DGeometry(
            self.minimum, nsteps=aabb_nsteps, step_vector=step_vector
        )
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

    @property
    def n_elements(self):
        return self.elements.shape[0]

    @property
    def n_nodes(self):
        return self.vertices.shape[0]

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

    def element_area(self, elements):
        tri_points = self.nodes[self.elements[elements, :], :]
        M_t = np.ones((tri_points.shape[0], 3, 3))
        M_t[:, :, 1:] = tri_points[:, :3, :]
        area = np.abs(np.linalg.det(M_t)) * 0.5
        return area

    def inside(self, pos):
        if pos.shape[1] > self.dimension:
            pos = pos[:, : self.dimension]

        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(self.dimension):
            inside *= pos[:, i] > self.minimum[None, i]
            inside *= pos[:, i] < self.maximum[None, i]
        return inside
