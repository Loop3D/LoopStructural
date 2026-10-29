"""The grid returned by ``feature.scalar_field`` must be registered to the values stored on it.

``scalar_field`` evaluates the feature at ``bounding_box.regular_grid()`` /
``bounding_box.cell_centres()`` and stores the results on a grid, so the grid's node
and cell positions must be the points the values were computed at. ``regular_grid`` puts
``nsteps`` points across the full extent (spacing ``node_step_vector``); ``step_vector`` is the
cell size used for interpolation supports and must not be used to describe that grid.

The feature is the analytical plane ``f(x, y, z) = x``, so the right value at any location is
its x coordinate and no interpolation is involved.
"""
import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.modelling.features._analytical_feature import (
    AnalyticalGeologicalFeature,
)


@pytest.fixture
def feature():
    # value == x coordinate
    return AnalyticalGeologicalFeature(
        name="x",
        vector=np.array([1.0, 0.0, 0.0]),
        origin=np.array([0.0, 0.0, 0.0]),
    )


@pytest.fixture
def bounding_box():
    return BoundingBox(origin=[0.0, 0.0, 0.0], maximum=[10.0, 10.0, 10.0], nsteps=[10, 10, 10])


def test_scalar_field_values_are_at_the_node_positions(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    nodes = np.asarray(grid.nodes)
    np.testing.assert_allclose(grid.properties["x"], nodes[:, 0], atol=1e-8)


def test_scalar_field_grid_spans_the_bounding_box(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    nodes = np.asarray(grid.nodes)
    np.testing.assert_allclose(nodes.min(axis=0), bounding_box.origin, atol=1e-8)
    np.testing.assert_allclose(nodes.max(axis=0), bounding_box.maximum, atol=1e-8)


def test_scalar_field_cell_values_are_at_the_cell_centres(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    centres = np.asarray(grid.cell_centres)
    np.testing.assert_allclose(grid.cell_properties["x"], centres[:, 0], atol=1e-8)


def test_cell_centres_are_inside_the_bounding_box(bounding_box):
    centres = bounding_box.cell_centres()
    assert np.all(centres >= bounding_box.origin - 1e-8)
    assert np.all(centres <= bounding_box.maximum + 1e-8)


@pytest.fixture
def offset_box():
    # non-zero origin and different node counts per axis: catches mixed-up axes and offsets
    return BoundingBox(origin=[100.0, -20.0, 5.0], maximum=[112.0, 4.0, 35.0], nsteps=[7, 5, 4])


def test_node_step_vector_is_the_spacing_of_regular_grid(offset_box):
    points = offset_box.regular_grid(local=False)
    for axis in range(3):
        coords = np.unique(np.round(points[:, axis], 10))
        assert len(coords) == offset_box.nsteps[axis]
        np.testing.assert_allclose(np.diff(coords), offset_box.node_step_vector[axis])


def test_cell_centres_count_and_midpoints(offset_box):
    centres = offset_box.cell_centres()
    assert centres.shape == (np.prod(offset_box.nsteps - 1), 3)
    for axis in range(3):
        coords = np.unique(np.round(centres[:, axis], 10))
        nodes = np.linspace(
            offset_box.origin[axis], offset_box.maximum[axis], offset_box.nsteps[axis]
        )
        np.testing.assert_allclose(coords, 0.5 * (nodes[:-1] + nodes[1:]))


def test_scalar_field_on_offset_anisotropic_box(feature, offset_box):
    feature = AnalyticalGeologicalFeature(
        name="f",
        vector=np.array([1.0, 2.0, 3.0]),
        origin=np.array([0.0, 0.0, 0.0]),
    )
    grid = feature.scalar_field(offset_box)
    expected_nodes = feature.evaluate_value(np.asarray(grid.nodes))
    np.testing.assert_allclose(grid.properties["f"], expected_nodes, atol=1e-8)
    expected_cells = feature.evaluate_value(np.asarray(grid.cell_centres))
    np.testing.assert_allclose(grid.cell_properties["f"], expected_cells, atol=1e-8)
