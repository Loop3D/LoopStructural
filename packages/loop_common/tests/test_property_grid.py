"""PropertyGrid: geometry comes from StructuredGrid3DGeometry, and cell lookups index cell data."""
import numpy as np
import pytest
from loop_common.geometry import PropertyGrid, StructuredGrid3DGeometry

GRIDS = {
    "unit_int": dict(origin=np.array([0, 0, 0]), step_vector=np.array([1, 1, 1]), nsteps=np.array([4, 4, 4])),
    "offset_anisotropic": dict(
        origin=np.array([10.0, -3.0, 2.0]),
        step_vector=np.array([0.5, 2.0, 1.5]),
        nsteps=np.array([7, 5, 4]),
    ),
    "minimal": dict(origin=np.zeros(3), step_vector=np.array([3.0, 2.0, 1.0]), nsteps=np.array([2, 2, 2])),
}


@pytest.fixture(params=list(GRIDS))
def grid(request):
    return PropertyGrid(**GRIDS[request.param])


def _axis_nodes(grid, axis):
    return grid.origin[axis] + grid.step_vector[axis] * np.arange(grid.nsteps[axis])


def test_geometry_is_built_from_the_current_fields(grid):
    assert isinstance(grid.geometry, StructuredGrid3DGeometry)
    np.testing.assert_array_equal(grid.geometry.nsteps, grid.nsteps)
    grid.nsteps = grid.nsteps + 1
    np.testing.assert_array_equal(grid.geometry.nsteps, grid.nsteps)


def test_maximum_and_counts(grid):
    np.testing.assert_allclose(grid.maximum, grid.origin + (grid.nsteps - 1) * grid.step_vector)
    np.testing.assert_array_equal(grid.nsteps_cells, grid.nsteps - 1)
    assert grid.n_nodes == np.prod(grid.nsteps)
    assert grid.n_cells == np.prod(grid.nsteps - 1)


def test_nodes_match_the_independent_construction(grid):
    xs, ys, zs = (_axis_nodes(grid, a) for a in range(3))
    expected = np.array(np.meshgrid(xs, ys, zs, indexing="ij")).reshape(3, -1, order="F").T
    np.testing.assert_allclose(grid.nodes, expected)


def test_cell_centres_match_the_independent_construction(grid):
    mids = [0.5 * (n[:-1] + n[1:]) for n in (_axis_nodes(grid, a) for a in range(3))]
    expected = np.array(np.meshgrid(*mids, indexing="ij")).reshape(3, -1, order="F").T
    assert grid.cell_centres.shape == (grid.n_cells, 3)
    np.testing.assert_allclose(grid.cell_centres, expected)


def test_cell_centres_are_inside_the_grid(grid):
    assert grid.inside(grid.cell_centres).all()
    assert (grid.cell_centres >= grid.origin).all() and (grid.cell_centres <= grid.maximum).all()


def test_each_cell_centre_is_in_its_own_cell(grid):
    np.testing.assert_array_equal(grid.cell_index_at(grid.cell_centres), np.arange(grid.n_cells))


def test_cell_value_at_reads_cell_data_in_storage_order(grid):
    # cell data defined as the x coordinate of the cell centre
    grid.cell_properties["x"] = grid.cell_centres[:, 0].copy()
    rng = np.random.default_rng(0)
    points = rng.uniform(grid.origin, grid.maximum, (200, 3))
    index = grid.cell_index_at(points)
    assert (index >= 0).all()
    # a point lies within half a step of the centre of the cell it is reported to be in
    centres = grid.cell_centres[index]
    assert (np.abs(points - centres) <= 0.5 * grid.step_vector + 1e-9).all()
    np.testing.assert_allclose(grid.cell_value_at("x", points), centres[:, 0])


def test_points_outside_give_minus_one_and_nan(grid):
    outside = np.array([grid.maximum + grid.step_vector, grid.origin - grid.step_vector])
    np.testing.assert_array_equal(grid.cell_index_at(outside), [-1, -1])
    grid.cell_properties["x"] = np.arange(grid.n_cells, dtype=float)
    assert np.isnan(grid.cell_value_at("x", outside)).all()
    assert not grid.inside(outside).any()


def test_boundary_points_belong_to_the_last_cell(grid):
    np.testing.assert_array_equal(grid.cell_index_at(grid.maximum), [grid.n_cells - 1])
    np.testing.assert_array_equal(grid.cell_index_at(grid.origin), [0])


def test_vector_valued_cell_data(grid):
    grid.cell_properties["v"] = np.column_stack([np.arange(grid.n_cells), -np.arange(grid.n_cells)]).astype(float)
    values = grid.cell_value_at("v", grid.cell_centres[[0, -1]])
    np.testing.assert_allclose(values, [[0, 0], [grid.n_cells - 1, -(grid.n_cells - 1)]])


# --------------------------------------------------------------------------- trilinear
LINEAR = np.array([1.5, -2.0, 0.75])


def _linear(points):
    return points @ LINEAR + 3.0


def test_trilinear_reproduces_a_linear_field_exactly(grid):
    grid.properties["f"] = _linear(grid.nodes)
    points = np.random.default_rng(1).uniform(grid.origin, grid.maximum, (300, 3))
    np.testing.assert_allclose(grid.evaluate_value("f", points), _linear(points), atol=1e-9)


def test_trilinear_returns_the_node_values_at_the_nodes(grid):
    values = np.random.default_rng(2).normal(size=grid.n_nodes)
    grid.properties["f"] = values
    np.testing.assert_allclose(grid.evaluate_value("f", grid.nodes), values, atol=1e-12)


def test_trilinear_is_nan_outside(grid):
    grid.properties["f"] = _linear(grid.nodes)
    outside = np.array([grid.maximum + grid.step_vector, grid.origin - grid.step_vector])
    assert np.isnan(grid.evaluate_value("f", outside)).all()


def test_trilinear_at_the_cell_midpoint_is_the_mean_of_its_corner_values(grid):
    values = np.random.default_rng(3).normal(size=grid.n_nodes)
    grid.properties["f"] = values
    geometry = grid.geometry
    centres = grid.cell_centres
    corners, _ = geometry.position_to_cell_corners(centres)
    np.testing.assert_allclose(grid.evaluate_value("f", centres), values[corners].mean(axis=1))


def test_trilinear_multicomponent(grid):
    grid.properties["v"] = np.column_stack([_linear(grid.nodes), -_linear(grid.nodes)])
    points = np.random.default_rng(4).uniform(grid.origin, grid.maximum, (50, 3))
    result = grid.evaluate_value("v", points)
    assert result.shape == (50, 2)
    np.testing.assert_allclose(result[:, 0], _linear(points), atol=1e-9)
    np.testing.assert_allclose(result[:, 1], -_linear(points), atol=1e-9)


def test_trilinear_agrees_with_the_interpolation_support(grid):
    from loop_common.supports import StructuredGrid

    support = StructuredGrid(
        origin=grid.origin, nsteps_cells=grid.nsteps - 1, step_vector=grid.step_vector
    )
    values = np.random.default_rng(5).normal(size=grid.n_nodes)
    grid.properties["f"] = values
    points = np.random.default_rng(6).uniform(grid.origin - grid.step_vector, grid.maximum + grid.step_vector, (200, 3))
    np.testing.assert_allclose(
        grid.evaluate_value("f", points), support.evaluate_value(points, values), equal_nan=True
    )


def test_transfer_a_property_onto_a_finer_support(grid):
    """The use case: evaluate a property at the nodes of another support."""
    from loop_common.supports import StructuredGrid

    grid.properties["f"] = _linear(grid.nodes)
    finer = StructuredGrid(
        origin=grid.origin,
        nsteps_cells=(grid.nsteps - 1) * 3,
        step_vector=grid.step_vector / 3.0,
    )
    transferred = grid.evaluate_value("f", finer.nodes)
    np.testing.assert_allclose(transferred, _linear(finer.nodes), atol=1e-9)
