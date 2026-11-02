"""Tests for the StructuredGrid lookup used by evaluate_value and evaluate_gradient.

The reference functions in this file are a frozen copy of the arithmetic that
StructuredGrid used before the common lookup was added (plan D6). They use only
numpy and the public geometry of the grid. They do not call _lookup_points or
_gradient_weights.

Float comparisons use rtol=1e-13 and atol=1e-14. The new code must do the same
operations in the same order, so the results must be equal to rounding level.
The small tolerance only allows for a different but equal expression.
Corner indexes and the inside mask are integer and boolean data and must be
exactly equal.
"""

import numpy as np
import pytest
from loop_common.supports import StructuredGrid

RTOL = 1e-13
ATOL = 1e-14

# Corner order of a cell: x changes fastest, then y, then z.
_XC = np.array([0, 1, 0, 1, 0, 1, 0, 1])
_YC = np.array([0, 0, 1, 1, 0, 0, 1, 1])
_ZC = np.array([0, 0, 0, 0, 1, 1, 1, 1])


# --------------------------------------------------------------------------
# Frozen reference calculations
# --------------------------------------------------------------------------
def _ref_inside(grid, pos):
    extent = np.maximum(np.abs(grid.maximum - grid.origin), 1.0)
    tol = np.finfo(float).eps * 32.0 * extent
    return np.all((pos >= (grid.origin - tol)) & (pos <= (grid.maximum + tol)), axis=1)


def _ref_cell_index(grid, pos):
    inside = _ref_inside(grid, pos)
    cell = np.zeros((pos.shape[0], 3), dtype=int)
    for axis in range(3):
        rel = pos[inside, axis] - grid.origin[axis]
        idx = np.floor(rel / grid.step_vector[axis]).astype(int)
        cell[inside, axis] = np.clip(idx, 0, grid.nsteps_cells[axis] - 1)
    return cell, inside


def _ref_corners(grid, pos):
    cell, inside = _ref_cell_index(grid, pos)
    nx, ny = grid.nsteps[0], grid.nsteps[1]
    base = cell[:, 0] + nx * cell[:, 1] + nx * ny * cell[:, 2]
    offsets = _XC + nx * _YC + nx * ny * _ZC
    corners = base[:, None] + offsets[None, :]
    corners[~inside] = -1
    return corners, inside


def _ref_local(grid, pos):
    cell, _ = _ref_cell_index(grid, pos)
    origin = grid.origin[None, :] + grid.step_vector[None, :] * cell
    return np.clip((pos - origin) / grid.step_vector[None, :], 0.0, 1.0)


def _ref_trilinear(lc):
    x, y, z = lc[:, 0], lc[:, 1], lc[:, 2]
    w = np.zeros((lc.shape[0], 8))
    w[:, 0] = (1 - x) * (1 - y) * (1 - z)
    w[:, 1] = x * (1 - y) * (1 - z)
    w[:, 2] = (1 - x) * y * (1 - z)
    w[:, 4] = (1 - x) * (1 - y) * z
    w[:, 5] = x * (1 - y) * z
    w[:, 6] = (1 - x) * y * z
    w[:, 3] = x * y * (1 - z)
    w[:, 7] = x * y * z
    return w


def _ref_gradient_weights(grid, lc):
    x, y, z = lc[:, 0], lc[:, 1], lc[:, 2]
    t = np.zeros((lc.shape[0], 3, 8))
    t[:, 0, 0] = (1 - z) * (y - 1)
    t[:, 0, 1] = (1 - y) * (1 - z)
    t[:, 0, 2] = -y * (1 - z)
    t[:, 0, 4] = -(1 - y) * z
    t[:, 0, 5] = (1 - y) * z
    t[:, 0, 6] = -y * z
    t[:, 0, 3] = y * (1 - z)
    t[:, 0, 7] = y * z
    t[:, 1, 0] = (x - 1) * (1 - z)
    t[:, 1, 1] = -x * (1 - z)
    t[:, 1, 2] = (1 - x) * (1 - z)
    t[:, 1, 4] = -(1 - x) * z
    t[:, 1, 5] = -x * z
    t[:, 1, 6] = (1 - x) * z
    t[:, 1, 3] = x * (1 - z)
    t[:, 1, 7] = x * z
    t[:, 2, 0] = -(1 - x) * (1 - y)
    t[:, 2, 1] = -x * (1 - y)
    t[:, 2, 2] = -(1 - x) * y
    t[:, 2, 4] = (1 - x) * (1 - y)
    t[:, 2, 5] = x * (1 - y)
    t[:, 2, 6] = (1 - x) * y
    t[:, 2, 3] = -x * y
    t[:, 2, 7] = x * y
    t[:, 0, :] /= grid.step_vector[0]
    t[:, 1, :] /= grid.step_vector[1]
    t[:, 2, :] /= grid.step_vector[2]
    return t


def _ref_value(grid, pos, prop):
    corners, inside = _ref_corners(grid, pos)
    v = np.full(corners.shape, np.nan)
    v[inside, :] = _ref_trilinear(_ref_local(grid, pos[inside, :]))
    v[inside, :] *= prop[corners[inside, :]]
    return np.sum(v, axis=1)


def _ref_gradient(grid, pos, prop):
    corners, inside = _ref_corners(grid, pos)
    t = np.zeros((corners.shape[0], 3, 8))
    t[inside] = _ref_gradient_weights(grid, _ref_local(grid, pos[inside, :]))
    for axis in range(3):
        t[inside, axis, :] *= prop[corners[inside, :]]
    return np.array(
        [np.sum(t[:, 0, :], axis=1), np.sum(t[:, 1, :], axis=1), np.sum(t[:, 2, :], axis=1)]
    ).T


# --------------------------------------------------------------------------
# Grids and points
# --------------------------------------------------------------------------
GRIDS = {
    "unit": dict(origin=[0.0, 0.0, 0.0], nsteps=[4, 5, 6], step_vector=[1.0, 1.0, 1.0]),
    "non_unit": dict(origin=[0.0, 0.0, 0.0], nsteps=[5, 3, 4], step_vector=[2.0, 0.5, 3.0]),
    "offset": dict(origin=[-10.5, 200.25, 7.0], nsteps=[3, 6, 4], step_vector=[1.5, 0.75, 2.5]),
}


def _make_grid(name):
    return StructuredGrid(**{k: np.array(v) for k, v in GRIDS[name].items()})


def _node_values(grid, seed=3):
    return np.random.default_rng(seed).normal(size=int(grid.n_nodes))


def _points(grid, seed=11):
    rng = np.random.default_rng(seed)
    lo, hi = grid.origin, grid.maximum
    span = hi - lo
    random_in = lo + rng.random((40, 3)) * span
    # Corners of the domain and of an interior cell.
    corners = np.array(
        [[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    )
    cell_nodes = lo + grid.step_vector * np.array([1.0, 1.0, 2.0])
    on_nodes = np.vstack([corners, cell_nodes[None, :]])
    # Points on faces and edges (one or two coordinates at a bound).
    faces = random_in[:6].copy()
    faces[0:2, 0] = lo[0]
    faces[2:4, 1] = hi[1]
    faces[4:6, 2] = hi[2]
    edges = random_in[6:10].copy()
    edges[0, [0, 1]] = [hi[0], lo[1]]
    edges[1, [1, 2]] = [hi[1], hi[2]]
    edges[2, [0, 2]] = [lo[0], hi[2]]
    edges[3, [0, 1]] = [hi[0], hi[1]]
    # Interior cell faces (points on a shared face between two cells).
    shared = random_in[10:14].copy()
    shared[:, 0] = lo[0] + grid.step_vector[0]
    eps = 1e-6 * span
    just_out = np.vstack([lo - eps, hi + eps, [lo[0] - eps[0], *(lo[1:] + 0.5 * span[1:])]])
    far_out = np.vstack([lo - 100.0 * span, hi + 100.0 * span])
    nan_rows = random_in[:3].copy()
    nan_rows[0, 0] = np.nan
    nan_rows[1, :] = np.nan
    nan_rows[2, 2] = np.nan
    return {
        "random": random_in,
        "nodes": on_nodes,
        "faces": faces,
        "edges": edges,
        "shared_face": shared,
        "just_outside": just_out,
        "far_outside": far_out,
        "nan": nan_rows,
        "mixed": np.vstack([random_in[:5], just_out, nan_rows, on_nodes[:3], far_out]),
        "empty": np.zeros((0, 3)),
    }


GRID_NAMES = list(GRIDS)
POINT_NAMES = list(_points(_make_grid("unit")))


@pytest.fixture(params=GRID_NAMES)
def grid(request):
    return _make_grid(request.param)


@pytest.fixture(params=POINT_NAMES)
def pts(request, grid):
    return _points(grid)[request.param]


# --------------------------------------------------------------------------
# Private helper contract (fail until _lookup_points / _gradient_weights exist)
# --------------------------------------------------------------------------
def test_lookup_points_matches_reference(grid, pts):
    corners, inside, local = grid._lookup_points(pts)
    ref_corners, ref_inside = _ref_corners(grid, pts)
    assert corners.shape == (pts.shape[0], 8)
    assert inside.shape == (pts.shape[0],)
    assert local.shape == (pts.shape[0], 3)
    assert inside.dtype == bool
    assert corners.dtype.kind == "i"
    assert local.dtype == np.float64
    np.testing.assert_array_equal(corners, ref_corners)
    np.testing.assert_array_equal(inside, ref_inside)
    np.testing.assert_allclose(local, _ref_local(grid, pts), rtol=RTOL, atol=ATOL, equal_nan=True)


def test_lookup_points_outside_rows(grid):
    pts = _points(grid)["mixed"]
    corners, inside, _local = grid._lookup_points(pts)
    assert not inside.all() and inside.any()
    # Outside points keep the corner marker -1. NaN points are outside.
    assert np.all(corners[~inside] == -1)
    assert np.all(corners[inside] >= 0)
    assert not inside[np.isnan(pts).any(axis=1)].any()


def test_lookup_points_corner_order_and_clipping(grid):
    # A point on the maximum corner belongs to the last cell, with local = 1.
    top = grid.maximum[None, :].copy()
    corners, inside, local = grid._lookup_points(top)
    assert inside[0]
    nx, ny = grid.nsteps[0], grid.nsteps[1]
    assert corners[0, 7] == grid.n_nodes - 1
    assert corners[0, 0] == corners[0, 7] - (1 + nx + nx * ny)
    np.testing.assert_array_equal(local, np.ones((1, 3)))
    # The corner order is x first, then y, then z.
    np.testing.assert_array_equal(corners[0] - corners[0, 0], _XC + nx * _YC + nx * ny * _ZC)


def test_lookup_points_accepts_single_point_and_list(grid):
    pt = grid.origin + 0.3 * (grid.maximum - grid.origin)
    c1, i1, l1 = grid._lookup_points(pt)
    c2, i2, l2 = grid._lookup_points([pt.tolist()])
    assert c1.shape == (1, 8)
    np.testing.assert_array_equal(c1, c2)
    np.testing.assert_array_equal(i1, i2)
    np.testing.assert_array_equal(l1, l2)


def test_gradient_weights_match_reference(grid):
    rng = np.random.default_rng(5)
    local = np.vstack(
        [
            rng.random((30, 3)),
            [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 1.0, 0.5], [np.nan, 0.5, 0.5]],
        ]
    )
    weights = grid._gradient_weights(local)
    assert isinstance(weights, np.ndarray)
    assert weights.shape == (local.shape[0], 3, 8)
    assert weights.dtype == np.float64
    # The weights include the division by the step vector (scale factors).
    np.testing.assert_allclose(
        weights, _ref_gradient_weights(grid, local), rtol=RTOL, atol=ATOL, equal_nan=True
    )


def test_gradient_weights_empty(grid):
    weights = grid._gradient_weights(np.zeros((0, 3)))
    assert weights.shape == (0, 3, 8)


def test_gradient_weights_match_element_gradient(grid):
    pts = _points(grid)["random"]
    local = grid._lookup_points(pts)[2]
    weights = grid._gradient_weights(local)
    old = grid.get_element_gradient_for_location(pts)[1]
    np.testing.assert_allclose(weights, old, rtol=RTOL, atol=ATOL)


# --------------------------------------------------------------------------
# Public methods against the reference (these freeze the present behaviour)
# --------------------------------------------------------------------------
def test_evaluate_value_matches_reference(grid, pts):
    prop = _node_values(grid)
    result = grid.evaluate_value(pts, prop)
    assert result.shape == (pts.shape[0],)
    assert result.dtype == np.float64
    expected = _ref_value(grid, pts, prop)
    np.testing.assert_allclose(result, expected, rtol=RTOL, atol=ATOL, equal_nan=True)
    # Outside points and NaN points give NaN.
    outside = ~_ref_inside(grid, pts)
    assert np.all(np.isnan(result[outside]))


def test_evaluate_gradient_matches_reference(grid, pts):
    prop = _node_values(grid)
    result = grid.evaluate_gradient(pts, prop)
    assert result.shape == (pts.shape[0], 3)
    assert result.dtype == np.float64
    expected = _ref_gradient(grid, pts, prop)
    np.testing.assert_allclose(result, expected, rtol=RTOL, atol=ATOL, equal_nan=True)
    # Outside points give zero gradient.
    outside = ~_ref_inside(grid, pts)
    assert np.all(result[outside] == 0.0)


def test_linear_field_value_and_gradient(grid):
    # A linear field is reproduced exactly by trilinear interpolation.
    nodes = grid.nodes
    coef = np.array([1.5, -2.0, 0.25])
    prop = nodes @ coef
    pts = _points(grid)["random"]
    np.testing.assert_allclose(grid.evaluate_value(pts, prop), pts @ coef, rtol=1e-10, atol=1e-9)
    np.testing.assert_allclose(
        grid.evaluate_gradient(pts, prop), np.tile(coef, (pts.shape[0], 1)), rtol=1e-9, atol=1e-9
    )


def test_property_array_size_error(grid):
    with pytest.raises(ValueError):
        grid.evaluate_value(np.zeros((1, 3)), np.zeros(3))
    with pytest.raises(ValueError):
        grid.evaluate_gradient(np.zeros((1, 3)), np.zeros(3))


# --------------------------------------------------------------------------
# Clipping of local coordinates inside the tolerance band (plan D6)
# --------------------------------------------------------------------------
# On a grid with a large extent the inside test accepts points up to
# tol = eps * 32 * extent beyond the box. For extent 1000 this is about 7e-12.
# A point 0.5 * tol beyond the box has an unclipped local coordinate that is
# outside [0, 1] by more than ATOL. Only the clip to [0, 1] gives exact 0 and 1.
CLIP_GRIDS = {
    "long_xy": dict(origin=[0.0, 0.0, 0.0], nsteps=[1000, 1000, 2], step_vector=[1.0, 1.0, 1.0]),
    "long_yz_offset": dict(
        origin=[-3.0, 500.5, 20.0], nsteps=[2, 1000, 1000], step_vector=[1.0, 1.0, 1.0]
    ),
}


def _tol(grid):
    # Same formula as the inside test of the grid.
    extent = np.maximum(np.abs(grid.maximum - grid.origin), 1.0)
    return np.finfo(float).eps * 32.0 * extent


def _band_points(grid):
    """Points in the tolerance band: (points, expected local coordinates)."""
    lo, hi, half = grid.origin, grid.maximum, 0.5 * _tol(grid)
    mid = lo + 0.5 * (hi - lo)
    pts = [hi + half, lo - half, np.where([True, False, True], hi + half, lo - half)]
    expected = [np.ones(3), np.zeros(3), np.array([1.0, 0.0, 1.0])]
    # One axis in the band, the other axes inside the box.
    for axis in range(3):
        for side, value in ((hi, 1.0), (lo, 0.0)):
            p = mid.copy()
            p[axis] = side[axis] + (half[axis] if value == 1.0 else -half[axis])
            e = np.full(3, np.nan)
            e[axis] = value
            pts.append(p)
            expected.append(e)
    return np.array(pts), np.array(expected)


@pytest.fixture(params=list(CLIP_GRIDS))
def clip_grid(request):
    return StructuredGrid(**{k: np.array(v) for k, v in CLIP_GRIDS[request.param].items()})


def test_band_points_overshoot_is_larger_than_atol(clip_grid):
    # Guard for the test data: the clip must matter at these points.
    assert np.max(_tol(clip_grid)) * 0.5 > 100 * ATOL


def test_band_points_are_inside(clip_grid):
    pts, _ = _band_points(clip_grid)
    assert np.all(clip_grid.inside(pts))
    assert np.all(_ref_inside(clip_grid, pts))


def test_position_to_local_coordinates_clips_in_tolerance_band(clip_grid):
    pts, expected = _band_points(clip_grid)
    local = clip_grid.position_to_local_coordinates(pts)
    assert np.all((local >= 0.0) & (local <= 1.0))
    known = ~np.isnan(expected)
    np.testing.assert_array_equal(local[known], expected[known])
    np.testing.assert_array_equal(local, _ref_local(clip_grid, pts))


def test_lookup_points_clips_in_tolerance_band(clip_grid):
    pts, expected = _band_points(clip_grid)
    corners, inside, local = clip_grid._lookup_points(pts)
    assert np.all(inside)
    assert np.all(corners >= 0)
    assert np.all((local >= 0.0) & (local <= 1.0))
    known = ~np.isnan(expected)
    np.testing.assert_array_equal(local[known], expected[known])
    np.testing.assert_array_equal(local, _ref_local(clip_grid, pts))


def test_evaluate_in_tolerance_band_matches_reference(clip_grid):
    pts, _ = _band_points(clip_grid)
    prop = _node_values(clip_grid)
    value = clip_grid.evaluate_value(pts, prop)
    assert not np.isnan(value).any()
    np.testing.assert_allclose(value, _ref_value(clip_grid, pts, prop), rtol=RTOL, atol=ATOL)
    gradient = clip_grid.evaluate_gradient(pts, prop)
    np.testing.assert_allclose(gradient, _ref_gradient(clip_grid, pts, prop), rtol=RTOL, atol=ATOL)


# --------------------------------------------------------------------------
# One cell-index lookup per call
# --------------------------------------------------------------------------
@pytest.fixture
def cell_index_counter(monkeypatch):
    calls = []
    original = StructuredGrid.position_to_cell_index

    def counting(self, pos):
        calls.append(np.shape(pos))
        return original(self, pos)

    monkeypatch.setattr(StructuredGrid, "position_to_cell_index", counting)
    return calls


@pytest.mark.parametrize("name", ["evaluate_value", "evaluate_gradient"])
def test_one_cell_index_lookup_per_call(name, cell_index_counter):
    grid = _make_grid("non_unit")
    prop = _node_values(grid)
    pts = _points(grid)["mixed"]
    getattr(grid, name)(pts, prop)
    assert len(cell_index_counter) == 1


@pytest.mark.parametrize("name", ["evaluate_value", "evaluate_gradient"])
def test_one_cell_index_lookup_for_inside_only_points(name, cell_index_counter):
    grid = _make_grid("offset")
    prop = _node_values(grid)
    pts = _points(grid)["random"]
    getattr(grid, name)(pts, prop)
    assert len(cell_index_counter) == 1


def test_lookup_points_uses_one_cell_index_lookup(cell_index_counter):
    grid = _make_grid("non_unit")
    grid._lookup_points(_points(grid)["mixed"])
    assert len(cell_index_counter) == 1
