"""
Tests for RectilinearGrid support.
"""

import numpy as np
import pytest
from loop_common.supports import RectilinearGrid

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def uniform_grid():
    """RectilinearGrid with uniform spacing — should behave like StructuredGrid."""
    x = np.linspace(0.0, 5.0, 6)
    y = np.linspace(0.0, 3.0, 4)
    z = np.linspace(0.0, 2.0, 3)
    return RectilinearGrid(x, y, z)


@pytest.fixture
def nonuniform_grid():
    """RectilinearGrid with deliberately non-uniform spacing."""
    x = np.array([0.0, 0.2, 0.5, 1.0, 2.0, 4.0, 5.0])
    y = np.array([0.0, 0.3, 0.7, 1.5, 3.0])
    z = np.array([0.0, 0.25, 0.75, 2.0])
    return RectilinearGrid(x, y, z)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_create_uniform(uniform_grid):
    assert uniform_grid is not None
    assert uniform_grid.n_nodes == 6 * 4 * 3
    assert uniform_grid.n_elements == 5 * 3 * 2


def test_create_nonuniform(nonuniform_grid):
    assert nonuniform_grid is not None
    assert nonuniform_grid.n_nodes == 7 * 5 * 4
    assert nonuniform_grid.n_elements == 6 * 4 * 3


def test_origin_maximum(nonuniform_grid):
    grid = nonuniform_grid
    assert np.allclose(grid.origin, [0.0, 0.0, 0.0])
    assert np.allclose(grid.maximum, [5.0, 3.0, 2.0])


def test_1d_arrays_required():
    with pytest.raises(ValueError, match="1-D"):
        RectilinearGrid(np.ones((3, 2)), np.linspace(0, 1, 3), np.linspace(0, 1, 3))


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------


def test_nodes_shape(nonuniform_grid):
    grid = nonuniform_grid
    assert grid.nodes.shape == (grid.n_nodes, 3)


def test_nodes_values(nonuniform_grid):
    grid = nonuniform_grid
    # Every xnodes value should appear in the x-column of nodes
    for xv in grid.xnodes:
        assert np.any(np.isclose(grid.nodes[:, 0], xv))
    for yv in grid.ynodes:
        assert np.any(np.isclose(grid.nodes[:, 1], yv))
    for zv in grid.znodes:
        assert np.any(np.isclose(grid.nodes[:, 2], zv))


# ---------------------------------------------------------------------------
# inside / position_to_cell_index
# ---------------------------------------------------------------------------


def test_barycentre_inside(nonuniform_grid):
    grid = nonuniform_grid
    assert np.all(grid.inside(grid.barycentre))


def test_outside_points_not_inside(nonuniform_grid):
    grid = nonuniform_grid
    outside_pts = np.array(
        [
            [-1.0, 1.0, 1.0],
            [6.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 4.0, 1.0],
        ]
    )
    assert not np.any(grid.inside(outside_pts))


def test_cell_index_range(nonuniform_grid):
    grid = nonuniform_grid
    pts = grid.barycentre
    idx, inside = grid.position_to_cell_index(pts)
    assert np.all(inside)
    assert np.all(idx[:, 0] < grid.nsteps_cells[0])
    assert np.all(idx[:, 1] < grid.nsteps_cells[1])
    assert np.all(idx[:, 2] < grid.nsteps_cells[2])
    assert np.all(idx >= 0)


def test_cell_index_correctness():
    """Each barycentre of cell (i,j,k) must land in cell (i,j,k)."""
    x = np.array([0.0, 1.0, 3.0, 6.0])
    y = np.array([0.0, 2.0, 5.0])
    z = np.array([0.0, 1.5, 4.0])
    grid = RectilinearGrid(x, y, z)
    centres = grid.barycentre
    idx, inside = grid.position_to_cell_index(centres)
    assert np.all(inside)
    gi_from_idx = grid.global_cell_indices(idx)
    gi_expected = np.arange(grid.n_elements)
    assert np.array_equal(gi_from_idx, gi_expected)


# ---------------------------------------------------------------------------
# Local coordinates
# ---------------------------------------------------------------------------


def test_local_coords_at_node_corners():
    """At node positions, local coords must be exactly 0 or 1."""
    x = np.array([0.0, 1.0, 3.0])
    y = np.array([0.0, 2.0])
    z = np.array([0.0, 0.5, 1.5])
    grid = RectilinearGrid(x, y, z)
    # Lower-left-front corner of the first cell → local = (0, 0, 0)
    pt_low = np.array([[0.0, 0.0, 0.0]])
    lc_low = grid.position_to_local_coordinates(pt_low)
    assert np.allclose(lc_low, 0.0)
    # Upper-right-back corner of the *last* cell → local = (1, 1, 1)
    pt_high = np.array([[3.0, 2.0, 1.5]])
    lc_high = grid.position_to_local_coordinates(pt_high)
    assert np.allclose(lc_high, 1.0)


def test_local_coords_midpoint():
    x = np.array([0.0, 2.0, 6.0])  # second cell has width 4
    y = np.array([0.0, 1.0])
    z = np.array([0.0, 1.0])
    grid = RectilinearGrid(x, y, z)
    # midpoint of second x-cell (x in [2,6]) at x=4 should give local_x=0.5
    pt = np.array([[4.0, 0.5, 0.5]])
    lc = grid.position_to_local_coordinates(pt)
    assert np.isclose(lc[0, 0], 0.5)


# ---------------------------------------------------------------------------
# DOF coefficients (trilinear partition of unity)
# ---------------------------------------------------------------------------


def test_dof_coefs_sum_to_one(nonuniform_grid):
    grid = nonuniform_grid
    pts = grid.barycentre
    coefs = grid.position_to_dof_coefs(pts)
    assert np.allclose(coefs.sum(axis=1), 1.0)


def test_dof_coefs_non_negative(nonuniform_grid):
    grid = nonuniform_grid
    coefs = grid.position_to_dof_coefs(grid.barycentre)
    assert np.all(coefs >= -1e-12)


# ---------------------------------------------------------------------------
# evaluate_value — interpolate a linear scalar field exactly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_evaluate_value_linear(nonuniform_grid, axis):
    """Trilinear interpolation must reproduce a linear field f = x_axis exactly."""
    grid = nonuniform_grid
    node_values = grid.nodes[:, axis]
    recovered = grid.evaluate_value(grid.barycentre, node_values)
    expected = grid.barycentre[:, axis]
    assert np.allclose(recovered, expected, atol=1e-10)


# ---------------------------------------------------------------------------
# evaluate_gradient — gradient of a linear field must be a unit basis vector
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_evaluate_gradient_linear(nonuniform_grid, axis):
    """Gradient of f = x_axis should be the axis unit vector everywhere."""
    grid = nonuniform_grid
    node_values = grid.nodes[:, axis]
    grads = grid.evaluate_gradient(grid.barycentre, node_values)
    expected = np.zeros((grid.n_elements, 3))
    expected[:, axis] = 1.0
    assert np.allclose(grads, expected, atol=1e-10)


# ---------------------------------------------------------------------------
# cell_centres
# ---------------------------------------------------------------------------


def test_cell_centres_inside(nonuniform_grid):
    grid = nonuniform_grid
    centres = grid.cell_centres(np.arange(grid.n_elements))
    assert np.all(grid.inside(centres))


def test_cell_centres_x_values():
    """Each cell centre x must equal the midpoint of that cell's x-interval."""
    x = np.array([0.0, 1.0, 3.0, 6.0])
    y = np.array([0.0, 1.0])
    z = np.array([0.0, 1.0])
    grid = RectilinearGrid(x, y, z)
    centres = grid.cell_centres(np.arange(grid.n_elements))
    np.tile([0.5, 2.0, 4.5], grid.n_elements // 3)
    assert np.allclose(np.sort(np.unique(centres[:, 0])), [0.5, 2.0, 4.5])


# ---------------------------------------------------------------------------
# build_scaled_operator_rows
# ---------------------------------------------------------------------------


def test_pure_second_derivative_shape(nonuniform_grid):
    grid = nonuniform_grid
    for axis in range(3):
        A, col, row = grid.build_scaled_operator_rows(axis)
        assert A.shape[1] == 3
        assert col.shape[1] == 3
        assert A.shape[0] == col.shape[0] == row.shape[0]


def test_mixed_second_derivative_shape(nonuniform_grid):
    grid = nonuniform_grid
    for ax, cx in [(0, 1), (0, 2), (1, 2)]:
        A, col, row = grid.build_scaled_operator_rows(ax, cx)
        assert A.shape[1] == 4
        assert col.shape[1] == 4
        assert A.shape[0] == col.shape[0] == row.shape[0]


def test_pure_operator_rows_sum_to_zero(nonuniform_grid):
    """Second-derivative stencil coefficients must sum to zero (constant field → zero curvature)."""
    grid = nonuniform_grid
    for axis in range(3):
        A, _, _ = grid.build_scaled_operator_rows(axis)
        assert np.allclose(A.sum(axis=1), 0.0, atol=1e-12)


def test_mixed_operator_rows_sum_to_zero(nonuniform_grid):
    grid = nonuniform_grid
    for ax, cx in [(0, 1), (0, 2), (1, 2)]:
        A, _, _ = grid.build_scaled_operator_rows(ax, cx)
        assert np.allclose(A.sum(axis=1), 0.0, atol=1e-12)


def test_pure_operator_uniform_matches_expected():
    """On a uniform grid, the pure d²/dx² stencil should give 1/h² [-1, 2, -1] (normalised)."""
    h = 0.5
    x = np.array([0.0, h, 2 * h, 3 * h])
    y = np.array([0.0, h])
    z = np.array([0.0, h])
    grid = RectilinearGrid(x, y, z)
    A, _, _ = grid.build_scaled_operator_rows(0)
    # Each row should be [2/(h*(2h)), -(4/h^2)/2, 2/(h*(2h))] = [1/h², -2/h², 1/h²]
    expected_coef = np.array([1 / h**2, -2 / h**2, 1 / h**2])
    assert np.allclose(A, expected_coef[None, :], atol=1e-10)


def test_global_indices_in_range(nonuniform_grid):
    grid = nonuniform_grid
    for axis in range(3):
        _, col, row = grid.build_scaled_operator_rows(axis)
        assert np.all(col >= 0)
        assert np.all(col < grid.n_nodes)
        assert np.all(row >= 0)
        assert np.all(row < grid.n_nodes)


# ---------------------------------------------------------------------------
# get_operators sentinel
# ---------------------------------------------------------------------------


def test_get_operators_returns_none_masks(nonuniform_grid):
    weights = dict.fromkeys(["dxx", "dyy", "dzz", "dxy", "dyz", "dxz"], 1.0)
    ops = nonuniform_grid.get_operators(weights)
    assert set(ops.keys()) == {"dxx", "dyy", "dzz", "dxy", "dyz", "dxz"}
    for name, (mask, _) in ops.items():
        assert mask is None, f"Expected None mask for operator '{name}'"


def test_get_operators_weight_values(nonuniform_grid):
    weights = {"dxx": 2.0, "dyy": 3.0, "dzz": 4.0, "dxy": 1.0, "dyz": 1.0, "dxz": 1.0}
    ops = nonuniform_grid.get_operators(weights)
    assert ops["dxx"][1] == 2.0
    assert ops["dyy"][1] == 3.0
    assert ops["dxy"][1] == pytest.approx(0.25)
