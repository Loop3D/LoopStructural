"""StructuredGrid2DGeometry: the geometry and node-value interpolation the 2D support builds on."""
import logging

import numpy as np
import pytest
from loop_common.geometry import StructuredGrid2DGeometry
from loop_common.supports import StructuredGrid2D

GRIDS = {
    "unit": dict(origin=[0.0, 0.0], nsteps=[5, 5], step_vector=[1.0, 1.0]),
    "aniso": dict(origin=[10.0, -3.0], nsteps=[7, 4], step_vector=[0.5, 2.0]),
}


@pytest.fixture(params=list(GRIDS))
def geometry(request):
    return StructuredGrid2DGeometry(**GRIDS[request.param])


def test_the_support_is_a_geometry():
    assert issubclass(StructuredGrid2D, StructuredGrid2DGeometry)


def test_elements_are_flat_node_indexes(geometry):
    assert geometry.elements.shape == (geometry.n_elements, 4)
    nodes = geometry.nodes
    # every element's four corners are the corners of one cell
    corners = nodes[geometry.elements]
    spans = corners.max(axis=1) - corners.min(axis=1)
    np.testing.assert_allclose(spans, np.broadcast_to(geometry.step_vector, spans.shape))


def test_index_round_trip(geometry):
    cells = np.arange(geometry.n_elements)
    ij = geometry.global_index_to_cell_index(cells)
    np.testing.assert_array_equal(geometry.global_cell_indices(ij), cells)
    nodes = np.arange(geometry.n_nodes)
    ij = geometry.global_index_to_node_index(nodes)
    np.testing.assert_array_equal(geometry.global_node_indices(ij), nodes)


def test_index_helpers_reject_bad_shapes(geometry):
    with pytest.raises(ValueError):
        geometry.global_cell_indices(np.array([0, 0]))
    with pytest.raises(ValueError):
        geometry.global_cell_indices(np.zeros((3, 3), dtype=int))


def test_bilinear_reproduces_a_linear_field(geometry):
    coefficients = np.array([1.5, -2.0])
    values = geometry.nodes @ coefficients + 3.0
    low = geometry.origin
    high = geometry.origin + (geometry.nsteps - 1) * geometry.step_vector
    points = np.random.default_rng(0).uniform(low + 1e-6, high - 1e-6, (200, 2))
    np.testing.assert_allclose(
        geometry.evaluate_value(points, values), points @ coefficients + 3.0, atol=1e-9
    )


def test_bilinear_returns_node_values_at_interior_nodes(geometry):
    values = np.random.default_rng(1).normal(size=geometry.n_nodes)
    nodes = geometry.nodes
    interior = (
        (nodes > geometry.origin + 1e-9).all(axis=1)
        & (nodes < geometry.origin + (geometry.nsteps - 1) * geometry.step_vector - 1e-9).all(axis=1)
    )
    np.testing.assert_allclose(
        geometry.evaluate_value(nodes[interior], values)[...], values[interior], atol=1e-12
    )


def test_print_geometry_works_on_a_2d_grid(geometry, caplog):
    with caplog.at_level(logging.INFO):
        geometry.print_geometry()
