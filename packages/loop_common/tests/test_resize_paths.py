"""The resize paths of the structured grids keep the extent and the counts.

Each resize setter must keep the grid extent, or enclose it, and must keep the node counts whole.
Setting the same value again and again must change nothing. The tests give a wrong answer on a
grid that drifts.
"""

import numpy as np
import pytest
from loop_common.geometry import StructuredGrid2DGeometry, StructuredGrid3DGeometry
from loop_common.supports import RectilinearGrid

# (origin, nodes per axis, step vector). nsteps counts nodes.
GRIDS_3D = {
    "G1": ((0.0, 0.0, 0.0), (11, 6, 4), (1.0, 2.0, 3.0)),
    "G2": ((100.0, -20.0, 5.0), (5, 5, 5), (0.5, 0.5, 0.5)),
    "G3": ((-1e4, 2e3, -30.0), (21, 3, 8), (10.0, 100.0, 7.5)),
}
GRIDS_2D = {
    "H1": ((0.0, 0.0), (5, 5), (1.0, 1.0)),
    "H2": ((10.0, -5.0), (8, 4), (0.5, 2.0)),
}
SETTERS = ("origin", "maximum", "nsteps", "step_vector")


def _grid_3d(name):
    origin, nsteps, step = GRIDS_3D[name]
    return StructuredGrid3DGeometry(origin=origin, nsteps=nsteps, step_vector=step)


def _state(grid):
    """Return a copy of the values that define the grid."""
    return {
        "origin": np.array(grid.origin, dtype=float),
        "maximum": np.array(grid.maximum, dtype=float),
        "step_vector": np.array(grid.step_vector, dtype=float),
        "nsteps": np.array(grid.nsteps),
    }


def _assert_state_equal(grid, expected):
    now = _state(grid)
    for key, value in expected.items():
        np.testing.assert_array_equal(now[key], value, err_msg=key)


def _tolerance(grid_state):
    """Length tolerance: 1e-12 * max(1, L) for each axis."""
    length = grid_state["maximum"] - grid_state["origin"]
    return 1e-12 * np.maximum(1.0, length)


@pytest.mark.parametrize("setter", SETTERS)
@pytest.mark.parametrize("name", list(GRIDS_3D))
def test_reapplying_same_value_10_times_changes_nothing_3d(name, setter):
    grid = _grid_3d(name)
    start = _state(grid)
    for _ in range(10):
        setattr(grid, setter, np.array(start[setter]))
    now = _state(grid)
    tol = _tolerance(start)
    for key in ("origin", "maximum", "step_vector"):
        assert np.all(np.abs(now[key] - start[key]) <= tol), (name, setter, key, now[key])
    np.testing.assert_array_equal(now["nsteps"], start["nsteps"])


@pytest.mark.parametrize("name", list(GRIDS_3D))
def test_nsteps_setter_keeps_extent(name):
    grid = _grid_3d(name)
    start = _state(grid)
    tol = _tolerance(start)
    for request in (np.array(start["nsteps"]) * 2, np.array(start["nsteps"])):
        grid.nsteps = request
        assert np.all(np.abs(grid.maximum - start["maximum"]) <= tol), (name, request, grid.maximum)
        assert grid.n_nodes == np.prod(request), (name, request, grid.n_nodes)


@pytest.mark.parametrize("name", list(GRIDS_3D))
def test_step_vector_setter_encloses_extent(name):
    grid = _grid_3d(name)
    start = _state(grid)
    new_step = start["step_vector"] * 1.5
    grid.step_vector = new_step.copy()
    tol = _tolerance(start)
    assert np.all(grid.maximum >= start["maximum"] - tol), (name, grid.maximum)
    assert np.all(grid.maximum < start["maximum"] + new_step), (name, grid.maximum)


@pytest.mark.parametrize("name", list(GRIDS_3D))
def test_step_vector_a_b_a_extent_drift(name):
    grid = _grid_3d(name)
    start = _state(grid)
    step_a = start["step_vector"].copy()
    step_b = step_a * 1.5
    grid.step_vector = step_a.copy()
    grid.step_vector = step_b.copy()
    grid.step_vector = step_a.copy()
    drift = grid.maximum - start["maximum"]
    print(f"{name}: A-B-A extent drift {drift}, limit {step_a + step_b}")
    tol = _tolerance(start)
    assert np.all(drift >= -tol), (name, drift)
    assert np.all(drift < step_a + step_b), (name, drift)


def test_nsteps_below_two_is_rejected_without_change():
    grid = _grid_3d("G1")
    start = _state(grid)
    with pytest.raises(ValueError):
        grid.nsteps = np.array([1, 6, 4])
    _assert_state_equal(grid, start)


def test_nsteps_not_whole_is_rejected_without_change():
    grid = _grid_3d("G1")
    start = _state(grid)
    for bad in ([5.5, 5, 5], [np.nan, 5, 5]):
        with pytest.raises(ValueError):
            grid.nsteps = np.array(bad)
        _assert_state_equal(grid, start)
    grid.nsteps = np.array([5.0, 5.0, 5.0])
    assert np.asarray(grid.nsteps).dtype.kind == "i"
    np.testing.assert_array_equal(grid.nsteps, [5, 5, 5])


@pytest.mark.parametrize("name", list(GRIDS_3D))
def test_set_nelements_keeps_extent(name):
    grid = _grid_3d(name)
    start = _state(grid)
    length = start["maximum"] - start["origin"]
    h = (np.prod(length) / 200) ** (1 / 3)
    result = grid.set_nelements(200)
    tol = _tolerance(start)
    assert np.all(np.abs(grid.maximum - start["maximum"]) <= tol), (name, grid.maximum)
    assert result == grid.n_elements
    expected = np.maximum(np.ceil(length / h).astype(int), 2)
    np.testing.assert_array_equal(grid.nsteps, expected)


def test_set_nelements_small_request_is_clamped():
    grid = _grid_3d("G1")
    start = _state(grid)
    grid.set_nelements(1)
    assert np.all(np.asarray(grid.nsteps) >= 2), grid.nsteps
    tol = _tolerance(start)
    assert np.all(np.abs(grid.maximum - start["maximum"]) <= tol), grid.maximum


@pytest.mark.parametrize("check", ("cells", "n_cell", "n_elements", "dtype", "maximum"))
@pytest.mark.parametrize("name", list(GRIDS_2D))
def test_2d_counts_are_derived(name, check):
    origin, nsteps, step = GRIDS_2D[name]
    grid = StructuredGrid2DGeometry(origin=origin, nsteps=nsteps, step_vector=step)
    grid.nsteps = [nsteps[0] + 3.0, nsteps[1] + 2.0]
    new = np.array([nsteps[0] + 3, nsteps[1] + 2])
    if check == "cells":
        np.testing.assert_array_equal(grid.nsteps_cells, new - 1)
    elif check == "n_cell":
        assert (grid.n_cell_x, grid.n_cell_y) == (new[0] - 1, new[1] - 1)
    elif check == "n_elements":
        assert grid.n_elements == (new[0] - 1) * (new[1] - 1)
    elif check == "dtype":
        assert np.asarray(grid.nsteps).dtype.kind == "i"
    else:
        # The 2D maximum keeps its legacy value: origin + nsteps * step_vector.
        np.testing.assert_allclose(grid.maximum, np.array(origin) + new * np.array(step))


@pytest.mark.parametrize("setter", ("nsteps", "step_vector", "origin"))
def test_rectilinear_resize_is_rejected(setter):
    grid = RectilinearGrid(np.linspace(0, 5, 6), np.linspace(0, 3, 4), np.linspace(0, 2, 3))
    nodes = np.array(grid.nodes)
    origin = np.array(grid.origin)
    nsteps = np.array(grid.nsteps)
    value = {
        "nsteps": np.array([7, 5, 4]),
        "step_vector": np.array([2.0, 2.0, 2.0]),
        "origin": np.array([1.0, 1.0, 1.0]),
    }[setter]
    with pytest.raises(NotImplementedError):
        setattr(grid, setter, value)
    np.testing.assert_array_equal(grid.nodes, nodes)
    np.testing.assert_array_equal(grid.origin, origin)
    np.testing.assert_array_equal(grid.nsteps, nsteps)
