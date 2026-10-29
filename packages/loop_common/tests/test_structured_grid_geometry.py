"""StructuredGrid3DGeometry must agree with the interpolation support it duplicates.

``loop_common.supports.StructuredGrid`` (via ``BaseStructuredSupport``) and
``StructuredGrid3DGeometry`` implement the same grid indexing. The tests below pin the shared
behaviour so the two cannot drift apart again; they compare the geometry class against the
support on the same inputs, and spell out the boundary rules.
"""
import numpy as np
import pytest
from loop_common.geometry import StructuredGrid3DGeometry
from loop_common.supports import StructuredGrid

# (origin, nodes per axis, step vector)
GRIDS = {
    "unit": ([0.0, 0.0, 0.0], [5, 5, 5], [1.0, 1.0, 1.0]),
    "offset_anisotropic": ([10.0, -3.0, 2.0], [7, 5, 4], [0.5, 2.0, 1.5]),
    "minimal": ([0.0, 0.0, 0.0], [2, 2, 2], [3.0, 2.0, 1.0]),
}


@pytest.fixture(params=list(GRIDS))
def pair(request):
    origin, nsteps, step = GRIDS[request.param]
    geometry = StructuredGrid3DGeometry(origin=origin, nsteps=nsteps, step_vector=step)
    support = StructuredGrid(origin=origin, nsteps_cells=np.array(nsteps) - 1, step_vector=step)
    return geometry, support


def _points(geometry, n=300, seed=0):
    rng = np.random.default_rng(seed)
    low = geometry.origin - 1.0
    high = geometry.maximum + 1.0
    return rng.uniform(low, high, (n, 3))


def test_counts_and_extent_agree(pair):
    geometry, support = pair
    np.testing.assert_array_equal(geometry.nsteps, support.nsteps)
    np.testing.assert_array_equal(geometry.nsteps_cells, support.nsteps_cells)
    np.testing.assert_allclose(geometry.maximum, support.maximum)
    assert geometry.n_nodes == support.n_nodes
    assert geometry.n_elements == support.n_elements
    np.testing.assert_allclose(geometry.nodes, support.nodes)


def test_position_to_cell_index_agrees(pair):
    geometry, support = pair
    points = _points(geometry)
    g_index, g_inside = geometry.position_to_cell_index(points)
    s_index, s_inside = support.position_to_cell_index(points)
    np.testing.assert_array_equal(g_inside, s_inside)
    np.testing.assert_array_equal(g_index, s_index)


def test_inside_agrees_including_the_boundary(pair):
    geometry, support = pair
    corners = np.array([geometry.origin, geometry.maximum])
    points = np.vstack([_points(geometry), corners])
    np.testing.assert_array_equal(geometry.inside(points), support.inside(points))
    assert geometry.inside(corners).all()


def test_index_conversions_agree(pair):
    geometry, support = pair
    cells = np.arange(geometry.n_elements)
    nodes = np.arange(geometry.n_nodes)
    np.testing.assert_array_equal(
        geometry.global_index_to_cell_index(cells), support.global_index_to_cell_index(cells)
    )
    np.testing.assert_array_equal(
        geometry.global_index_to_node_index(nodes), support.global_index_to_node_index(nodes)
    )
    cell_ijk = geometry.global_index_to_cell_index(cells)
    node_ijk = geometry.global_index_to_node_index(nodes)
    np.testing.assert_array_equal(
        geometry.global_cell_indices(cell_ijk), support.global_cell_indices(cell_ijk)
    )
    np.testing.assert_array_equal(
        geometry.global_node_indices(node_ijk), support.global_node_indices(node_ijk)
    )


def test_elements_are_flat_node_indexes(pair):
    geometry, support = pair
    assert geometry.elements.shape == (geometry.n_elements, 8)
    np.testing.assert_array_equal(geometry.elements, support.elements)


def test_cell_index_round_trips_through_the_flat_index(pair):
    geometry, _ = pair
    cells = np.arange(geometry.n_elements)
    ijk = geometry.global_index_to_cell_index(cells)
    np.testing.assert_array_equal(geometry.global_cell_indices(ijk), cells)


def test_point_on_the_maximum_boundary_is_in_the_last_cell():
    geometry = StructuredGrid3DGeometry(origin=[0, 0, 0], nsteps=[5, 5, 5], step_vector=[1, 1, 1])
    index, inside = geometry.position_to_cell_index(np.array([[4.0, 4.0, 4.0], [0.0, 0.0, 0.0]]))
    assert inside.all()
    np.testing.assert_array_equal(index, [[3, 3, 3], [0, 0, 0]])


def test_point_outside_is_flagged_not_indexed():
    geometry = StructuredGrid3DGeometry(origin=[0, 0, 0], nsteps=[5, 5, 5], step_vector=[1, 1, 1])
    index, inside = geometry.position_to_cell_index(np.array([[6.0, 2.0, 2.0], [-0.5, 2.0, 2.0]]))
    assert not inside.any()
    np.testing.assert_array_equal(index, [[0, 0, 0], [0, 0, 0]])


def test_single_point_is_accepted():
    geometry = StructuredGrid3DGeometry(origin=[0, 0, 0], nsteps=[5, 5, 5], step_vector=[1, 1, 1])
    index, inside = geometry.position_to_cell_index(np.array([1.5, 2.5, 0.5]))
    np.testing.assert_array_equal(index, [[1, 2, 0]])
    assert inside.tolist() == [True]


def test_index_helpers_reject_bad_shapes():
    geometry = StructuredGrid3DGeometry(origin=[0, 0, 0], nsteps=[3, 3, 3], step_vector=[1, 1, 1])
    with pytest.raises(ValueError):
        geometry.global_cell_indices(np.array([0, 0, 0]))
    with pytest.raises(ValueError):
        geometry.global_cell_indices(np.zeros((4, 2), dtype=int))
