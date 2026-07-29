"""
P2TetMesh based on cartesian grid for piecewise quadratic interpolation.

P2 tetrahedra have 10 nodes:
- 4 corner nodes (vertices of the tetrahedron)
- 6 edge midpoint nodes (one at the center of each edge)

This mesh adds mid-edge nodes to a structured cartesian grid.
"""

import numpy as np
from ._3d_base_structured import BaseStructuredSupport
from . import SupportType
from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


class P2TetMesh(BaseStructuredSupport):
    """P2 (piecewise quadratic) tetrahedral mesh from a structured cartesian grid.
    
    P2 elements have 10 nodes per tetrahedron (4 vertices + 6 edge midpoints).
    This class builds a mesh with both vertex and edge midpoint nodes.
    """

    def __init__(self, origin=np.zeros(3), nsteps=np.ones(3) * 10, step_vector=np.ones(3)):
        BaseStructuredSupport.__init__(self, origin, nsteps, step_vector)
        self.type = SupportType.P2StructuredTetMesh
        
        # P1 tetra masks (same as TetMesh - used for vertex connectivity)
        self.tetra_mask_even = np.array(
            [[7, 1, 2, 4], [6, 2, 4, 7], [5, 1, 4, 7], [0, 1, 2, 4], [3, 1, 2, 7]]
        )
        self.tetra_mask = np.array(
            [[0, 6, 5, 3], [7, 3, 5, 6], [4, 0, 5, 6], [2, 0, 3, 6], [1, 0, 3, 5]]
        )
        self.hessian = np.array(
            [
                [
                    [4, 4, 0, 0, 0, 0, -8, 0, 0, 0],
                    [4, 0, 0, 0, 0, 0, -4, 4, 0, -4],
                    [4, 0, 0, 0, 0, -4, -4, 0, 4, 0],
                ],
                [
                    [4, 0, 0, 0, 0, 0, -4, 4, 0, -4],
                    [4, 0, 4, 0, 0, 0, 0, 0, 0, -8],
                    [4, 0, 0, 0, 4, -4, 0, 0, 0, -4],
                ],
                [
                    [4, 0, 0, 0, 0, -4, -4, 0, 4, 0],
                    [4, 0, 0, 0, 4, -4, 0, 0, 0, -4],
                    [4, 0, 0, 4, 0, -8, 0, 0, 0, 0],
                ],
            ]
        )
        
        self.cg = None
        self._elements = None
        self._nodes = None
        self._edge_node_indices = None
        self._p1_mesh = None

    def onGeometryChange(self):
        """Reset cached properties when geometry changes."""
        self._elements = None
        self._nodes = None
        self._edge_node_indices = None
        self._p1_mesh = None
        if self.interpolator is not None:
            self.interpolator.reset()

    def _get_p1_mesh(self):
        """Return a structured P1 mesh with the same geometry."""
        if self._p1_mesh is None:
            from ._3d_structured_tetra import TetMesh

            self._p1_mesh = TetMesh(self.origin, self.nsteps - 1, self.step_vector)
        return self._p1_mesh

    @property
    def n_nodes(self) -> int:
        """Total number of nodes: vertices + edge midpoints."""
        return self.n_vertices + self.n_edge_nodes

    @property
    def ntetra(self) -> int:
        """Number of tetrahedra: 5 per cell."""
        return np.prod(self.nsteps_cells) * 5

    @property
    def n_elements(self) -> int:
        """Number of elements (same as ntetra)."""
        return self.ntetra

    @property
    def n_cells(self) -> int:
        """Number of cells in the grid."""
        return np.prod(self.nsteps_cells)

    @property
    def n_vertices(self) -> int:
        """Number of vertex nodes (corners of the grid)."""
        return np.prod(self.nsteps)

    @property
    def n_edge_nodes(self) -> int:
        """Number of edge midpoint nodes."""
        if self._edge_node_indices is None:
            self._build_edge_node_map()
        return len(self._edge_node_indices)

    @property
    def nodes(self):
        """All nodes: vertices followed by edge midpoints."""
        if self._nodes is None:
            self._nodes = self._compute_all_nodes()
        return self._nodes

    def _compute_all_nodes(self) -> np.ndarray:
        """Compute vertex and edge midpoint node coordinates."""
        x = np.arange(self.nsteps[0])
        y = np.arange(self.nsteps[1])
        z = np.arange(self.nsteps[2])
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        vertex_indexes = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        vertex_nodes = self.node_indexes_to_position(vertex_indexes)
        
        # Compute edge midpoint nodes
        edge_nodes = self._compute_edge_nodes()
        
        # Concatenate vertices and edge nodes
        all_nodes = np.vstack([vertex_nodes, edge_nodes])
        return all_nodes

    def _compute_edge_nodes(self) -> np.ndarray:
        """Compute coordinates of all edge midpoint nodes."""
        if self._edge_node_indices is None:
            self._build_edge_node_map()

        x = np.arange(self.nsteps[0])
        y = np.arange(self.nsteps[1])
        z = np.arange(self.nsteps[2])
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        vertex_indexes = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        vertex_nodes = self.node_indexes_to_position(vertex_indexes)
        edge_nodes_list = [None] * len(self._edge_node_indices)
        for (v1_idx, v2_idx), edge_index in self._edge_node_indices.items():
            pos1 = vertex_nodes[v1_idx]
            pos2 = vertex_nodes[v2_idx]
            edge_nodes_list[edge_index - self.n_vertices] = 0.5 * (pos1 + pos2)

        return np.array(edge_nodes_list)

    def _build_edge_node_map(self):
        """Build a mapping from vertex pairs to edge node indices."""
        if self._edge_node_indices is not None:
            return self._edge_node_indices
        
        edge_map = {}
        n_vertices = self.n_vertices

        from ._3d_structured_tetra import TetMesh
        dummy_mesh = TetMesh(self.origin, self.nsteps - 1, self.step_vector)
        p1_tetras = dummy_mesh.get_elements()

        for tetra in p1_tetras:
            tetra = np.asarray(tetra, dtype=int)
            tetra_edges = (
                (tetra[0], tetra[1]),
                (tetra[0], tetra[2]),
                (tetra[0], tetra[3]),
                (tetra[1], tetra[2]),
                (tetra[1], tetra[3]),
                (tetra[2], tetra[3]),
            )

            for edge in tetra_edges:
                key = tuple(sorted((int(edge[0]), int(edge[1]))))
                if key not in edge_map:
                    edge_map[key] = n_vertices + len(edge_map)
        
        self._edge_node_indices = edge_map
        return edge_map

    @property
    def elements(self):
        """Get all P2 elements (10-node tetrahedra)."""
        if self._elements is None:
            self._elements = self.get_elements()
        return self._elements

    def get_elements(self):
        """Get all P2 tetrahedra with 10 nodes each.
        
        Returns
        -------
        np.ndarray
            Array of shape (n_tetras, 10) with node indices for each P2 tetrahedron.
        """
        edge_map = self._build_edge_node_map()
        tetras_list = []
        
        # Get P1 tetrahedra (vertex connectivity)
        from ._3d_structured_tetra import TetMesh
        # TetMesh.__init__ adds 1 to nsteps, so pass nsteps - 1 here.
        dummy_mesh = TetMesh(self.origin, self.nsteps - 1, self.step_vector)
        p1_tetras = dummy_mesh.get_elements()
        # For each P1 tetrahedron, add the 6 edge midpoint nodes
        for p1_tetra_vertices in p1_tetras:
            # p1_tetra_vertices is a list of 4 vertex indices
            p2_tetra = list(p1_tetra_vertices)
            
            # Add the 6 edge midpoint nodes
            edges = [
                (p1_tetra_vertices[2], p1_tetra_vertices[3]),
                (p1_tetra_vertices[0], p1_tetra_vertices[3]),
                (p1_tetra_vertices[0], p1_tetra_vertices[1]),
                (p1_tetra_vertices[1], p1_tetra_vertices[2]),
                (p1_tetra_vertices[1], p1_tetra_vertices[3]),
                (p1_tetra_vertices[0], p1_tetra_vertices[2]),
            ]
            
            for v1, v2 in edges:
                edge_key = tuple(sorted([v1, v2]))
                if edge_key in edge_map:
                    p2_tetra.append(edge_map[edge_key])
                else:
                    # This shouldn't happen in a well-formed mesh
                    raise ValueError(f"Edge {edge_key} not found in edge map")
            
            tetras_list.append(p2_tetra)
        
        return np.array(tetras_list, dtype=np.int64)

    def evaluate_shape(self, locations: np.ndarray):
        """Evaluate quadratic tetrahedral shape functions at locations."""
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        N = np.zeros((c.shape[0], 10))

        for i in range(c.shape[1]):
            N[:, i] = (2 * c[:, i] - 1) * c[:, i]

        N[:, 4] = 4 * c[:, 3] * c[:, 2]
        N[:, 5] = 4 * c[:, 0] * c[:, 3]
        N[:, 6] = 4 * c[:, 0] * c[:, 1]
        N[:, 7] = 4 * c[:, 1] * c[:, 2]
        N[:, 8] = 4 * c[:, 1] * c[:, 3]
        N[:, 9] = 4 * c[:, 0] * c[:, 2]

        return N, elements, inside

    def evaluate_shape_derivatives(self, locations: np.ndarray, elements=None):
        """Evaluate quadratic tetrahedral shape derivatives at locations."""
        locations = np.array(locations)
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(locations)
        else:
            M = np.ones((elements.shape[0], 4, 4))
            M[:, :, 1:] = self.nodes[self.elements[elements], :][:, :4, :]
            points_ = np.ones((locations.shape[0], 4))
            points_[:, 1:] = locations
            minv = np.linalg.inv(M)
            c = np.einsum("lij,li->lj", minv, points_)
            verts = self.nodes[self.elements[elements][:, :4]]

        jac = np.array(
            [
                [
                    (verts[:, 1, 0] - verts[:, 0, 0]),
                    (verts[:, 1, 1] - verts[:, 0, 1]),
                    (verts[:, 1, 2] - verts[:, 0, 2]),
                ],
                [
                    (verts[:, 2, 0] - verts[:, 0, 0]),
                    (verts[:, 2, 1] - verts[:, 0, 1]),
                    (verts[:, 2, 2] - verts[:, 0, 2]),
                ],
                [
                    (verts[:, 3, 0] - verts[:, 0, 0]),
                    (verts[:, 3, 1] - verts[:, 0, 1]),
                    (verts[:, 3, 2] - verts[:, 0, 2]),
                ],
            ]
        )
        r = c[:, 1]
        s = c[:, 2]
        t = c[:, 3]
        jac = np.swapaxes(jac, 0, 2)
        dN = np.zeros((elements.shape[0], 3, 10))

        dN[:, 0, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 0, 1] = 4 * r - 1
        dN[:, 0, 2] = 0
        dN[:, 0, 3] = 0
        dN[:, 0, 4] = 0
        dN[:, 0, 5] = -4 * t
        dN[:, 0, 6] = -8 * r - 4 * s - 4 * t + 4
        dN[:, 0, 7] = 4 * s
        dN[:, 0, 8] = 4 * t
        dN[:, 0, 9] = -4 * s

        dN[:, 1, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 1, 1] = 0
        dN[:, 1, 2] = 4 * s - 1
        dN[:, 1, 3] = 0
        dN[:, 1, 4] = 4 * t
        dN[:, 1, 5] = -4 * t
        dN[:, 1, 6] = -4 * r
        dN[:, 1, 7] = 4 * r
        dN[:, 1, 8] = 0
        dN[:, 1, 9] = -4 * r - 8 * s - 4 * t + 4

        dN[:, 2, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 2, 1] = 0
        dN[:, 2, 2] = 0
        dN[:, 2, 3] = 4 * t - 1
        dN[:, 2, 4] = 4 * s
        dN[:, 2, 5] = -4 * r - 4 * s - 8 * t + 4
        dN[:, 2, 6] = -4 * r
        dN[:, 2, 7] = 0
        dN[:, 2, 8] = 4 * r
        dN[:, 2, 9] = -4 * s

        d_n = np.linalg.inv(jac)
        d_n = d_n.swapaxes(1, 2)
        d_n = d_n @ dN
        return d_n, elements

    def evaluate_shape_d2(self, indexes: np.ndarray) -> np.ndarray:
        """Evaluate second derivatives of the P2 tetrahedral shape functions."""
        vertices = self.nodes[self.elements[indexes], :]

        jac = np.array(
            [
                [
                    (vertices[:, 1, 0] - vertices[:, 0, 0]),
                    (vertices[:, 1, 1] - vertices[:, 0, 1]),
                    (vertices[:, 1, 2] - vertices[:, 0, 2]),
                ],
                [
                    (vertices[:, 2, 0] - vertices[:, 0, 0]),
                    (vertices[:, 2, 1] - vertices[:, 0, 1]),
                    (vertices[:, 2, 2] - vertices[:, 0, 2]),
                ],
                [
                    (vertices[:, 3, 0] - vertices[:, 0, 0]),
                    (vertices[:, 3, 1] - vertices[:, 0, 1]),
                    (vertices[:, 3, 2] - vertices[:, 0, 2]),
                ],
            ]
        )
        jac = jac.swapaxes(0, 2)
        jac = jac.swapaxes(1, 2)
        jac = np.linalg.inv(jac)

        d2 = np.zeros((vertices.shape[0], 6, self.elements.shape[1]))
        ii = 0
        for i in range(3):
            for j in range(i, 3):
                for k in range(3):
                    for l in range(3):
                        d2[:, ii, :] += (
                            jac[:, i, k, None] * jac[:, j, l, None] * self.hessian[None, k, l, :]
                        )
                ii += 1
        return d2

    def get_quadrature_points(self, npts: int = 3):
        """Return face quadrature points for shared tetra faces."""
        if npts not in (1, 3):
            raise ValueError("Only 1-point and 3-point quadrature are supported")

        vertices = self.nodes[self.shared_elements]
        if npts == 3:
            cp = np.zeros((vertices.shape[0], 3, 3))
            reference_points = np.array([[1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])

            cp[:, 0, :] = (
                vertices[:, 0, :] * (1 - reference_points[0, 0] - reference_points[1, 0])
                + vertices[:, 1, :] * (reference_points[0, 0])
                + vertices[:, 2, :] * (reference_points[1, 0])
            )
            cp[:, 1, :] = (
                vertices[:, 0, :] * (1 - reference_points[0, 1] - reference_points[1, 1])
                + vertices[:, 1, :] * (reference_points[0, 1])
                + vertices[:, 2, :] * (reference_points[1, 1])
            )
            cp[:, 2, :] = (
                vertices[:, 0, :] * (1 - reference_points[0, 2] - reference_points[1, 2])
                + vertices[:, 1, :] * (reference_points[0, 2])
                + vertices[:, 2, :] * (reference_points[1, 2])
            )
            weights = np.zeros((vertices.shape[0], 3))
            weights[:, :] = 1 / 6
            return cp, weights

        cp = np.zeros((vertices.shape[0], 1, 3))
        reference_points = np.array([[1 / 3], [1 / 3]])
        cp[:, 0, :] = (
            vertices[:, 0, :] * (1 - reference_points[0, 0] - reference_points[1, 0])
            + vertices[:, 1, :] * (reference_points[0, 0])
            + vertices[:, 2, :] * (reference_points[1, 0])
        )
        weights = np.zeros((vertices.shape[0], 1))
        weights[:, :] = 1 / 2
        return cp, weights

    @property
    def shared_element_relationships(self):
        """Shared face relationships from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().shared_element_relationships

    @property
    def shared_elements(self):
        """Shared face node indices from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().shared_elements

    @property
    def shared_element_norm(self):
        """Face normals from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().shared_element_norm

    @property
    def shared_element_size(self):
        """Face areas from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().shared_element_size

    @property
    def shared_element_scale(self):
        """Scaled face areas from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().shared_element_scale

    @property
    def neighbours(self):
        """Neighbour relationships from the underlying P1 tetra mesh."""
        return self._get_p1_mesh().neighbours

    @property
    def element_size(self):
        """Calculate the volume of tetrahedra using the 4 corner vertices."""
        vecs = (
            self.nodes[self.elements[:, 1:4], :]
            - self.nodes[self.elements[:, 0, None], :]
        )
        return np.abs(np.linalg.det(vecs)) / 6

    @property
    def barycentre(self) -> np.ndarray:
        """Return the barycentres of all tetrahedra."""
        tetra_vertices = self.elements[:, :4]
        barycentre = np.sum(self.nodes[tetra_vertices][:, :, :], axis=1) / 4.0
        return barycentre

    def evaluate_value(self, pos: np.ndarray, property_array: np.ndarray) -> np.ndarray:
        """Evaluate value at given positions using P2 interpolation."""
        values = np.zeros(pos.shape[0])
        values[:] = np.nan

        N, tetras, inside = self.evaluate_shape(pos)
        values[inside] = np.sum(
            N[inside, :] * property_array[self.elements[tetras[inside], :]], axis=1
        )
        return values

    def evaluate_gradient(self, pos: np.ndarray, property_array: np.ndarray) -> np.ndarray:
        """Evaluate the gradient of an interpolant at the locations.
        
        This uses the quadratic tetrahedral basis on the 10 P2 nodes.
        """
        values = np.zeros(pos.shape)
        values[:] = np.nan

        element_gradients, tetras = self.evaluate_shape_derivatives(pos)
        _, _, inside = self.evaluate_shape(pos)
        values[inside, :] = np.einsum(
            "ijk,ik->ij",
            element_gradients[inside, :, :],
            property_array[self.elements[tetras[inside], :]],
        )

        return values

    def get_element_for_location(self, pos: np.ndarray):
        """Determine the tetrahedron from a numpy array of points.
        
        This uses only the vertex nodes for locating points.
        """
        pos = np.array(pos)
        pos = pos[:, : self.dimension]
        inside = self.inside(pos)
        
        # Create vertex-only elements for point location (P1 mesh)
        vertices = np.zeros((pos.shape[0], 5, 4, 3))
        vertices[:] = np.nan
        
        cell_indexes, inside = self.position_to_cell_index(pos)
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        
        corner_indexes = self.cell_corner_indexes(cell_indexes)
        vert_positions = self.node_indexes_to_position(corner_indexes)
        
        vertices[even_mask, :, :, :] = vert_positions[even_mask, :, :][:, self.tetra_mask_even, :]
        vertices[~even_mask, :, :, :] = vert_positions[~even_mask, :, :][:, self.tetra_mask, :]
        
        vap = pos[:, None, :] - vertices[:, :, 0, :]
        vbp = pos[:, None, :] - vertices[:, :, 1, :]
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
        
        mask = np.all(c >= 0, axis=2)
        i, j = np.where(mask)
        pairs = dict(zip(i, j))
        mask[:] = False
        mask[list(pairs.keys()), list(pairs.values())] = True
        
        inside = np.logical_and(inside, np.any(mask, axis=1))
        
        even_mask = np.sum(cell_indexes, axis=1) % 2 == 0
        gi = self.global_node_indices(corner_indexes)
        
        tetras = np.zeros((corner_indexes.shape[0], 5, 4)).astype(int)
        tetras[even_mask, :, :] = gi[even_mask, :][:, self.tetra_mask_even]
        tetras[~even_mask, :, :] = gi[~even_mask, :][:, self.tetra_mask]
        
        inside = np.logical_and(inside, self.inside(pos))
        
        vertices_return = np.zeros((pos.shape[0], 4, 3))
        vertices_return[:] = np.nan
        mask[~inside, :] = False
        vertices_return[inside, :, :] = vertices[mask, :, :]
        
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

    def get_element_gradient_for_location(self, pos: np.ndarray):
        """Get the gradient of the tetra for a location.
        
        Uses vertex-only evaluation (P1 gradients).
        """
        vertices, bc, tetras, inside = self.get_element_for_location(pos)
        ps = vertices
        
        m = np.array([
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
        ])
        
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.zeros_like(m)
        element_gradients[:] = np.nan
        element_gradients[inside, :, :] = np.linalg.inv(m[inside, :, :])
        
        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        
        return vertices, element_gradients, tetras, inside

    def tetra_global_index(self, indices, tetra_index):
        """Get the global index of a tetra from cell index and local tetra index."""
        return (
            tetra_index
            + indices[:, 0] * 5
            + self.nsteps_cells[0] * indices[:, 1] * 5
            + self.nsteps_cells[0] * self.nsteps_cells[1] * indices[:, 2] * 5
        )

    def inside(self, pos: np.ndarray):
        """Check if points are inside the mesh domain."""
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside
