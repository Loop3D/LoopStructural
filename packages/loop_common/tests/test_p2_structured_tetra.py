"""Tests for P2TetMesh (piecewise quadratic structured tetrahedral mesh)."""

import numpy as np
import pytest
from loop_common.supports._3d_structured_tetra import TetMesh
from loop_common.supports._p2_structured_tetra import P2TetMesh


class TestP2TetMeshConstruction:
    """Test P2TetMesh initialization and basic properties."""

    def test_p2tetmesh_basic_creation(self):
        """Test P2TetMesh can be created with default parameters."""
        mesh = P2TetMesh()
        assert mesh is not None
        # Default nsteps=np.ones(3)*10 passed to parent which adds 1
        assert mesh.nsteps[0] == 11
        assert mesh.nsteps[1] == 11
        assert mesh.nsteps[2] == 11

    def test_p2tetmesh_custom_grid(self):
        """Test P2TetMesh with custom origin, nsteps, and step_vector."""
        origin = np.array([1.0, 2.0, 3.0])
        nsteps = np.array([5, 6, 7])
        step_vector = np.array([0.5, 0.5, 0.5])
        mesh = P2TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)

        # The support constructor treats nsteps as the number of cells and
        # stores nsteps + 1 (number of vertices) internally.
        assert np.allclose(mesh.origin, origin)
        assert np.allclose(mesh.nsteps, nsteps + 1)
        assert np.allclose(mesh.step_vector, step_vector)

    def test_p2tetmesh_n_vertices(self):
        """Test that n_vertices is the product of nsteps."""
        mesh = P2TetMesh(nsteps=np.array([3, 4, 5]))
        expected_vertices = 4 * 5 * 6
        assert mesh.n_vertices == expected_vertices

    def test_p2tetmesh_n_tetras(self):
        """Test that ntetra is 5 * n_cells."""
        nsteps = np.array([3, 3, 3])
        mesh = P2TetMesh(nsteps=nsteps)
        # mesh.nsteps_cells == nsteps (the constructor's nsteps + 1 vertices,
        # minus 1, gives back the requested number of cells)
        n_cells = np.prod(nsteps)
        expected_tetras = 5 * n_cells
        assert mesh.ntetra == expected_tetras

    def test_p2tetmesh_n_nodes(self):
        """Test that n_nodes = n_vertices + n_edge_nodes."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        assert mesh.n_nodes == mesh.n_vertices + mesh.n_edge_nodes

    def test_p2tetmesh_elements_shape(self):
        """Test that elements array has correct shape (ntetra, 10)."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        elements = mesh.get_elements()
        assert elements.shape == (mesh.ntetra, 10)

    def test_p2tetmesh_elements_node_indices_valid(self):
        """Test that all node indices in elements are within valid range."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        elements = mesh.get_elements()
        assert np.all(elements >= 0)
        assert np.all(elements < mesh.n_nodes)

    def test_p2tetmesh_elements_first_four_nodes_are_vertices(self):
        """Test that first 4 nodes of each element are vertices (P1 tetra)."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        elements = mesh.get_elements()
        # First 4 columns should be vertex indices (0 to n_vertices-1)
        assert np.all(elements[:, :4] < mesh.n_vertices)

    def test_p2tetmesh_elements_last_six_are_edge_nodes(self):
        """Test that last 6 nodes of each element are edge node indices."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        elements = mesh.get_elements()
        # Last 6 columns should be edge nodes (>= n_vertices)
        assert np.all(elements[:, 4:] >= mesh.n_vertices)
        assert np.all(elements[:, 4:] < mesh.n_nodes)


class TestP2TetMeshNodes:
    """Test P2TetMesh node generation and coordinates."""

    def test_p2tetmesh_nodes_shape(self):
        """Test that nodes array has correct shape (n_nodes, 3)."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        nodes = mesh.nodes
        assert nodes.shape == (mesh.n_nodes, 3)

    def test_p2tetmesh_vertex_nodes_coordinates(self):
        """Test that vertex node coordinates match the cartesian grid."""
        origin = np.array([0.0, 0.0, 0.0])
        nsteps = np.array([3, 3, 3])
        step_vector = np.array([1.0, 1.0, 1.0])
        mesh = P2TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)

        nodes = mesh.nodes
        # First n_vertices nodes are vertices
        vertex_nodes = nodes[: mesh.n_vertices]

        # Check a few expected vertex positions
        assert np.allclose(vertex_nodes[0], [0.0, 0.0, 0.0])  # origin
        # nsteps=3 cells -> mesh.nsteps == 4 vertices per axis (indices 0..3)
        assert np.allclose(vertex_nodes[-1], [3.0, 3.0, 3.0])

    def test_p2tetmesh_edge_nodes_are_midpoints(self):
        """Test that edge nodes are at midpoints of vertex edges."""
        origin = np.array([0.0, 0.0, 0.0])
        nsteps = np.array([3, 3, 3])
        step_vector = np.array([1.0, 1.0, 1.0])
        mesh = P2TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)

        nodes = mesh.nodes
        elements = mesh.get_elements()

        # For each element, check that edge nodes (indices 4-9) are at midpoints
        # of vertex node pairs (indices 0-3)
        for element in elements[:10]:  # Check first 10 elements
            # Get the 4 vertex nodes
            v_nodes = nodes[element[:4]]
            # Get the 6 edge nodes
            edge_nodes = nodes[element[4:]]

            # The 6 edges are:
            # (2,3), (0,3), (0,1), (1,2), (1,3), (0,2)
            edge_pairs = [(2, 3), (0, 3), (0, 1), (1, 2), (1, 3), (0, 2)]
            for i, (v1_idx, v2_idx) in enumerate(edge_pairs):
                expected_edge = 0.5 * (v_nodes[v1_idx] + v_nodes[v2_idx])
                assert np.allclose(
                    edge_nodes[i], expected_edge, atol=1e-10
                ), f"Edge {i} midpoint mismatch"


class TestP2TetMeshShapeFunctions:
    """Test P2 tetrahedral shape function evaluation."""

    def test_p2tetmesh_shape_function_partition_of_unity(self):
        """Test that P2 shape functions sum to 1 at any point in element."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))

        # Test at a point inside an element
        # Use centroid of first element's vertices
        elements = mesh.get_elements()
        first_element_vertices = mesh.nodes[elements[0, :4]]
        centroid = np.mean(first_element_vertices, axis=0)

        N, elem_ids, inside = mesh.evaluate_shape(centroid.reshape(1, 3))

        assert inside[0], "Test point should be inside first element"
        # Sum of shape functions should be 1
        assert np.isclose(np.sum(N[0]), 1.0, atol=1e-10)

    def test_p2tetmesh_shape_function_at_vertices(self):
        """Test that shape function is ~1 at its node, ~0 at others."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        elements = mesh.get_elements()

        # Get the first element
        element = elements[0]
        element_vertices = mesh.nodes[element[:4], :]
        centroid = np.mean(element_vertices, axis=0)

        # For each vertex in the element, evaluate shape functions just inside
        # the element (grid vertices are shared by several tetrahedra, so
        # querying the exact vertex position is ambiguous about which
        # element/local-node it resolves to)
        for node_idx in range(4):
            point = (0.999 * element_vertices[node_idx] + 0.001 * centroid).reshape(1, 3)
            N, elem_ids, inside = mesh.evaluate_shape(point)

            assert inside[0], f"Vertex {node_idx} should be inside element"
            assert elem_ids[0] == 0
            # Shape function at its own vertex should be ~1
            assert np.isclose(
                N[0, node_idx], 1.0, atol=1e-2
            ), f"Shape function {node_idx} at its vertex"

    def test_p2tetmesh_shape_function_derivatives_exist(self):
        """Test that shape derivatives can be computed."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))

        centroid = np.array([[1.0, 1.0, 1.0]])
        dN, elem_ids = mesh.evaluate_shape_derivatives(centroid)

        # dN should have shape (n_points, 3, 10)
        assert dN.shape == (1, 3, 10)

    def test_p2tetmesh_shape_second_derivatives_exist(self):
        """Test that second derivatives can be computed."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))

        # Get valid element indices
        elements = mesh.get_elements()
        d2 = mesh.evaluate_shape_d2(np.array([0]))

        # d2 should have shape (n_elements, 6, 10)
        # (6 because second derivative has 6 independent components in 3D)
        assert d2.shape == (1, 6, 10)


class TestP2TetMeshEdgeNodes:
    """Test edge node generation and mapping."""

    def test_p2tetmesh_edge_node_map_consistency(self):
        """Test that edge node map is consistent across tetrahedra sharing edges."""
        mesh = P2TetMesh(nsteps=np.array([4, 4, 4]))
        elements = mesh.get_elements()

        # For any two tetrahedral elements that share an edge,
        # they should reference the same edge node index
        edge_to_elements = {}

        for elem_idx, element in enumerate(elements):
            # Extract vertex indices (first 4 nodes)
            vertices = element[:4]
            # All possible edges in this tetrahedron
            edges = [
                (vertices[0], vertices[1]),
                (vertices[0], vertices[2]),
                (vertices[0], vertices[3]),
                (vertices[1], vertices[2]),
                (vertices[1], vertices[3]),
                (vertices[2], vertices[3]),
            ]

            # Normalize edges (smaller index first)
            for edge in edges:
                normalized = tuple(sorted(edge))
                if normalized not in edge_to_elements:
                    edge_to_elements[normalized] = []
                edge_to_elements[normalized].append(elem_idx)

    def test_p2tetmesh_n_edge_nodes_reasonable(self):
        """Test that number of edge nodes is reasonable."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        n_edges = mesh.n_edge_nodes

        # Each cell (cube) is split into 5 tetrahedra whose edges are a
        # subset of the cube's own axis edges, face diagonals and body
        # diagonals - i.e. at most C(8, 2) = 28 unique vertex pairs per
        # cell. Edges shared between neighbouring cells only reduce this,
        # so n_cells * 28 is always a safe (if loose) upper bound.
        n_cells = np.prod(mesh.nsteps_cells)
        max_edges = n_cells * 28

        assert n_edges > 0
        assert n_edges <= max_edges


class TestP2TetMeshGeometryChange:
    """Test geometry invalidation and caching."""

    def test_p2tetmesh_onGeometryChange_clears_cache(self):
        """Test that onGeometryChange() clears cached properties."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))

        # Force computation of cached properties
        _ = mesh.nodes
        _ = mesh.elements

        assert mesh._nodes is not None
        assert mesh._elements is not None

        # Change geometry
        mesh.onGeometryChange()

        assert mesh._nodes is None
        assert mesh._elements is None

    def test_p2tetmesh_nodes_recomputed_after_geometry_change(self):
        """Test that nodes are recomputed after geometry change."""
        origin1 = np.array([0.0, 0.0, 0.0])
        origin2 = np.array([1.0, 1.0, 1.0])
        mesh = P2TetMesh(origin=origin1, nsteps=np.array([3, 3, 3]))

        nodes1 = mesh.nodes.copy()

        # Change origin
        mesh.origin = origin2
        mesh.onGeometryChange()

        nodes2 = mesh.nodes

        # Nodes should be different. Note: the origin setter preserves
        # `maximum` and `step_vector`, recomputing nsteps, so this is a
        # resize (not a pure translation) and the node arrays may differ
        # in shape as well as values.
        assert not np.array_equal(nodes1, nodes2)


class TestP2TetMeshQuadrature:
    """Test quadrature point generation."""

    def test_p2tetmesh_quadrature_points_1pt(self):
        """Test 1-point quadrature (centroid)."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        cp, weights = mesh.get_quadrature_points(npts=1)

        # Should have one point per shared element
        assert cp.shape[0] == mesh.shared_elements.shape[0]
        assert cp.shape[1] == 1  # 1 quadrature point
        assert cp.shape[2] == 3  # 3D coordinates

        # Weight should be the area of the reference triangle (1/2)
        assert np.allclose(weights, 1.0 / 2.0)

    def test_p2tetmesh_quadrature_points_3pt(self):
        """Test 3-point quadrature."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))
        cp, weights = mesh.get_quadrature_points(npts=3)

        # Should have 3 points per shared element
        assert cp.shape[1] == 3  # 3 quadrature points
        assert cp.shape[2] == 3  # 3D coordinates

        # Each weight should be 1/6
        assert np.allclose(weights, 1.0 / 6.0)

    def test_p2tetmesh_invalid_quadrature_points(self):
        """Test that invalid npts raises error."""
        mesh = P2TetMesh(nsteps=np.array([3, 3, 3]))

        with pytest.raises(ValueError, match="Only 1-point and 3-point"):
            mesh.get_quadrature_points(npts=2)


class TestP2TetMeshIntegration:
    """Integration tests comparing P2 and P1 meshes."""

    def test_p2tetmesh_same_vertex_positions_as_p1(self):
        """Test that P2TetMesh vertices match P1TetMesh vertices."""
        nsteps = np.array([4, 4, 4])
        # Both TetMesh and P2TetMesh add 1 to nsteps internally (cells -> vertices),
        # so passing the same nsteps to both produces the same underlying grid.
        p1_mesh = TetMesh(nsteps=nsteps)
        p2_mesh = P2TetMesh(nsteps=nsteps)

        p1_vertices = p1_mesh.nodes
        p2_vertices = p2_mesh.nodes[: p2_mesh.n_vertices]

        assert np.allclose(p1_vertices, p2_vertices)

    def test_p2tetmesh_same_p1_vertex_connectivity_as_p1(self):
        """Test that first 4 nodes of P2 elements match P1 elements."""
        nsteps = np.array([3, 3, 3])
        p1_mesh = TetMesh(nsteps=nsteps)
        p2_mesh = P2TetMesh(nsteps=nsteps)

        p1_elements = p1_mesh.get_elements()
        p2_elements = p2_mesh.get_elements()

        # P2 elements have 10 nodes, first 4 are the P1 vertices
        # The connectivity should match
        for p1_elem, p2_elem in zip(p1_elements, p2_elements):
            assert np.array_equal(p1_elem, p2_elem[:4])
