"""BoundingBox.structured_grid returns a PropertyGrid whose nodes are the regular_grid points."""
import numpy as np
import pytest
from loop_common.geometry import BoundingBox, PropertyGrid

BOXES = {
    "unit": dict(origin=[0.0, 0.0, 0.0], maximum=[1.0, 1.0, 1.0], nsteps=[10, 10, 10]),
    "offset_anisotropic": dict(origin=[100.0, -20.0, 5.0], maximum=[112.0, 4.0, 35.0], nsteps=[7, 5, 4]),
}


@pytest.fixture(params=list(BOXES))
def box(request):
    return BoundingBox(**BOXES[request.param])


def test_returns_a_property_grid(box):
    assert isinstance(box.structured_grid(), PropertyGrid)


def test_nodes_are_the_regular_grid_points(box):
    grid = box.structured_grid()
    np.testing.assert_allclose(grid.nodes, box.regular_grid(local=False, order="F"))
    np.testing.assert_array_equal(grid.nsteps, box.nsteps)


def test_grid_spans_the_bounding_box(box):
    grid = box.structured_grid()
    np.testing.assert_allclose(grid.origin, box.origin)
    np.testing.assert_allclose(grid.maximum, box.maximum)


def test_cell_centres_are_the_bounding_box_cell_centres(box):
    grid = box.structured_grid()
    np.testing.assert_allclose(grid.cell_centres, box.reproject(box.cell_centres(order="F")))


def test_data_is_stored_and_copied(box):
    nodes = np.arange(float(box.nsteps.prod()))
    cells = np.arange(float((box.nsteps - 1).prod()))
    node_data = {"n": nodes}
    grid = box.structured_grid(cell_data={"c": cells}, vertex_data=node_data, name="g")
    assert grid.name == "g"
    np.testing.assert_array_equal(grid.properties["n"], nodes)
    np.testing.assert_array_equal(grid.cell_properties["c"], cells)
    grid.properties["other"] = nodes
    assert "other" not in node_data
    other = box.structured_grid()
    assert other.properties == {} and other.cell_properties == {}


def test_values_computed_on_the_box_land_on_the_grid_nodes(box):
    """The property that matters for scalar_field: a value computed at regular_grid() is at that node."""
    grid = box.structured_grid()
    grid.properties["x"] = box.regular_grid(local=False, order="F")[:, 0]
    np.testing.assert_allclose(grid.properties["x"], grid.nodes[:, 0])
    np.testing.assert_allclose(grid.evaluate_value("x", grid.nodes), grid.nodes[:, 0])


def test_local_coordinates_grid_contains_the_rotated_box():
    box = BoundingBox(origin=[0.0, 0.0, 0.0], maximum=[10.0, 10.0, 10.0], nsteps=[5, 5, 5])
    angle = np.radians(45.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]]
    )
    box.set_local_transform(local_origin=[0.0, 0.0, 0.0], rotation_matrix=rotation)
    corners = box.project(box.corners)
    grid = box.structured_grid(local_coordinates=True)
    np.testing.assert_allclose(grid.origin, corners.min(axis=0), atol=1e-8)
    np.testing.assert_allclose(grid.maximum, corners.max(axis=0), atol=1e-8)
    np.testing.assert_array_equal(grid.nsteps, box.nsteps)
