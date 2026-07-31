"""Pure unstructured mesh geometry: nodes/elements/neighbours containers."""

import numpy as np
from scipy import sparse

from ._aabb import _initialise_aabb
from ._face_table import _init_face_table
from ._structured_grid_2d import StructuredGrid2DGeometry
from ._structured_grid_3d import StructuredGrid3DGeometry


class UnstructuredMeshGeometry:
    """An unstructured tetrahedral mesh defined by nodes, elements and neighbours."""

    dimension = 3

    def __init__(self, nodes: np.ndarray, elements: np.ndarray, neighbours: np.ndarray, aabb_nsteps=None):
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
            self.aabb_grid = StructuredGrid3DGeometry(self.minimum, nsteps=[2, 2, 2], step_vector=[1, 1, 1])
        else:
            if aabb_nsteps is None:
                box_vol = np.prod(self.maximum - self.minimum)
                element_volume = box_vol / (len(self.elements) / 20)
                step_vector = np.zeros(3)
                step_vector[:] = element_volume ** (1.0 / 3.0)
                aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(int)
                aabb_nsteps[aabb_nsteps < 2] = 2
            aabb_nsteps = np.array(aabb_nsteps, dtype=int)
            step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
            self.aabb_grid = StructuredGrid3DGeometry(self.minimum, nsteps=aabb_nsteps, step_vector=step_vector)
        self._aabb_table = sparse.csr_matrix((self.aabb_grid.n_elements, len(self.elements)), dtype=bool)
        self._shared_element_relationships = np.zeros((self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int)
        self._shared_elements = np.zeros((self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int)

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
        return self.neighbours

    @property
    def shared_element_norm(self):
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        v2 = self.nodes[elements[:, 2], :] - self.nodes[elements[:, 0], :]
        return np.cross(v1, v2, axisa=1, axisb=1)

    @property
    def shared_element_size(self):
        norm = self.shared_element_norm
        return 0.5 * np.linalg.norm(norm, axis=1)

    @property
    def element_size(self):
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
    """An unstructured triangular mesh defined by vertices, elements and neighbours."""

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
            step_vector = np.zeros(2)
            step_vector[:] = element_volume ** (1.0 / 2.0)
            aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(int)
            aabb_nsteps[aabb_nsteps < 2] = 2
        step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
        self.aabb_grid = StructuredGrid2DGeometry(self.minimum, nsteps=aabb_nsteps, step_vector=step_vector)
        self._aabb_table = sparse.csr_matrix((self.aabb_grid.n_elements, len(self.elements)), dtype=bool)
        self._shared_element_relationships = np.zeros((self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int)
        self._shared_elements = np.zeros((self.neighbours[self.neighbours >= 0].flatten().shape[0], self.dimension), dtype=int)

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
        return self.elements.shape[1]

    @property
    def nodes(self):
        return self.vertices

    @property
    def barycentre(self):
        element_idx = np.arange(0, self.n_elements)
        elements = self.elements[element_idx]
        barycentre = np.sum(self.nodes[elements][:, :3, :], axis=1) / 3.0
        return barycentre

    @property
    def shared_element_norm(self):
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        norm = np.zeros_like(v1)
        norm[:, 0] = v1[:, 1]
        norm[:, 1] = -v1[:, 0]
        return norm
