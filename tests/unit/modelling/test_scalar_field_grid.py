"""The grid returned by ``feature.scalar_field`` must be registered to the values stored on it.

``scalar_field`` evaluates the feature at ``bounding_box.regular_grid()`` /
``bounding_box.cell_centres()`` and stores the results on a grid built from
``bounding_box.origin``, ``step_vector`` and ``nsteps``. Those two descriptions of the box
disagree (``step_vector`` is ``extent / nsteps`` but ``regular_grid`` places ``nsteps`` points
across the full extent), so the stored values are not at the node positions the grid reports.

The feature is the analytical plane ``f(x, y, z) = x``, so the right value at any location is
its x coordinate and no interpolation is involved.

These are strict xfails: they fail today. When the bounding box geometry is made consistent
they will XPASS, which fails the run - remove the markers then.
"""
import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.modelling.features._analytical_feature import (
    AnalyticalGeologicalFeature,
)

GEOMETRY_BUG = pytest.mark.xfail(
    strict=True,
    reason="BoundingBox.step_vector (extent / nsteps) is inconsistent with regular_grid "
    "(nsteps points across the full extent), so scalar_field grids are misregistered",
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


@GEOMETRY_BUG
def test_scalar_field_values_are_at_the_node_positions(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    nodes = np.asarray(grid.nodes)
    np.testing.assert_allclose(grid.properties["x"], nodes[:, 0], atol=1e-8)


@GEOMETRY_BUG
def test_scalar_field_grid_spans_the_bounding_box(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    nodes = np.asarray(grid.nodes)
    np.testing.assert_allclose(nodes.min(axis=0), bounding_box.origin, atol=1e-8)
    np.testing.assert_allclose(nodes.max(axis=0), bounding_box.maximum, atol=1e-8)


@GEOMETRY_BUG
def test_scalar_field_cell_values_are_at_the_cell_centres(feature, bounding_box):
    grid = feature.scalar_field(bounding_box)
    centres = np.asarray(grid.cell_centres)
    np.testing.assert_allclose(grid.cell_properties["x"], centres[:, 0], atol=1e-8)


@GEOMETRY_BUG
def test_cell_centres_are_inside_the_bounding_box(bounding_box):
    centres = bounding_box.cell_centres()
    assert np.all(centres >= bounding_box.origin - 1e-8)
    assert np.all(centres <= bounding_box.maximum + 1e-8)
