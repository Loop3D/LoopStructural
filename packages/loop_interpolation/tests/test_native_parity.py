"""Parity tests: native solver against the Python matrix solver (T027).

Each test builds one problem twice with ``InterpolatorBuilder``. One build uses
``solver="cg"`` (matrix solver, tight tolerance). The other build uses
``solver="native"``. The two fields must differ by at most 1e-6 of the range.

All tests skip when the native core is not installed.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("loop_interpolation_core")

from loop_interpolation import InterpolatorBuilder, InterpolatorType, StructuredGrid  # noqa: E402

# Number of nodes on each axis (30 x 30 x 30 = 27000 nodes).
N_NODES = 30
SEED = 20240607
REL_TOL = 1e-6
CG_KWARGS = {"rtol": 1e-12, "atol": 1e-12, "maxiter": 200000}
# The native solver gets the same tolerance and iteration limit as the reference.
NATIVE_KWARGS = {"tol": 1e-12, "maxiter": 50000}

# Small background weight for the regularisation terms that are not under test.
# It keeps the problem well posed when one term is tested alone.
BACKGROUND = 0.05
ALL_REGULARISATION = ("dxx", "dyy", "dzz", "dxy", "dxz", "dyz")
ZERO_DATA_WEIGHTS = {"cpw": 0.0, "gpw": 0.0, "npw": 0.0, "tpw": 0.0, "ipw": 0.0}


def _make_grid(n_nodes: int = N_NODES) -> StructuredGrid:
    # nsteps_cells is the number of cells, so the node count is n_nodes.
    return StructuredGrid(
        origin=np.zeros(3),
        nsteps_cells=np.full(3, n_nodes - 1),
        step_vector=np.ones(3),
    )


def _field(p: np.ndarray) -> np.ndarray:
    x, y, z = (p[:, i] / (N_NODES - 1) for i in range(3))
    return 3.0 * x + 1.5 * np.sin(2.0 * y) + 2.0 * z * z + 0.5 * x * y


def _field_gradient(p: np.ndarray) -> np.ndarray:
    x, y, z = (p[:, i] / (N_NODES - 1) for i in range(3))
    g = np.column_stack([3.0 + 0.5 * y, 3.0 * np.cos(2.0 * y) + 0.5 * x, 4.0 * z])
    return g / (N_NODES - 1)


def _points(rng, n: int) -> np.ndarray:
    # All data is inside the grid, away from the border.
    return rng.uniform(2.0, N_NODES - 3.0, size=(n, 3))


def _data(family: str) -> dict:
    """Return the constraint arrays for one family (set-method name -> array)."""
    rng = np.random.default_rng(SEED)
    value_pts = _points(rng, 80)
    value = np.column_stack([value_pts, _field(value_pts), rng.uniform(0.5, 2.0, 80)])
    if family == "value":
        return {"add_value_constraints": value}
    if family == "gradient":
        p = _points(rng, 30)
        g = _field_gradient(p)
        g /= np.linalg.norm(g, axis=1)[:, None]
        return {"add_value_constraints": value, "add_gradient_constraints": np.column_stack([p, g])}
    if family == "normal":
        p = _points(rng, 30)
        n = _field_gradient(p)
        return {
            "add_value_constraints": value,
            "add_normal_constraints": np.column_stack([p, n, np.ones(len(p))]),
        }
    if family == "tangent":
        p = _points(rng, 30)
        t = np.cross(_field_gradient(p), np.tile([0.0, 0.0, 1.0], (len(p), 1)))
        t /= np.linalg.norm(t, axis=1)[:, None]
        return {"add_value_constraints": value, "add_tangent_constraints": np.column_stack([p, t])}
    if family == "interface":
        p = _points(rng, 12)
        ids = np.repeat(np.arange(4, dtype=float), 3)
        return {
            "add_value_constraints": value,
            "set_interface_constraints": np.column_stack([p, ids]),
        }
    # Regularisation and border families use value data only.
    return {"add_value_constraints": value}


def _build(solver: str, data: dict, setup_kwargs: dict) -> np.ndarray:
    builder = InterpolatorBuilder(
        interpolatortype=InterpolatorType.FINITE_DIFFERENCE,
        support=_make_grid(),
    )
    if solver == "native":
        builder.use_solver("native", **NATIVE_KWARGS)
    else:
        builder.use_solver(solver, **CG_KWARGS)
    for method, array in data.items():
        if hasattr(builder, method):
            getattr(builder, method)(array)
        else:
            # The builder has no adder for this family. Use the interpolator.
            getattr(builder.interpolator, method)(array)
    builder.setup_interpolator(interface_pair_mode="star", **setup_kwargs)
    builder.solve()
    interpolator = builder.build()
    assert interpolator.up_to_date
    return np.asarray(interpolator.c, dtype=float).copy()


def _assert_parity(data: dict, setup_kwargs: dict) -> None:
    reference = _build("cg", data, setup_kwargs)
    native = _build("native", data, setup_kwargs)
    assert native.shape == reference.shape
    span = float(np.ptp(reference))
    assert span > 0.0
    assert np.max(np.abs(native - reference)) <= REL_TOL * span


def _background_regularisation(**overrides) -> dict:
    kwargs = {name: BACKGROUND for name in ALL_REGULARISATION}
    kwargs.update(overrides)
    return kwargs


@pytest.mark.parametrize("family", ["value", "gradient", "normal", "tangent", "interface"])
def test_native_matches_cg_for_constraint_family(family):
    # Default regularisation, one data family (plus the value data that fixes the level).
    weights = {**ZERO_DATA_WEIGHTS, "cpw": 1.0}
    if family != "value":
        weights[{"gradient": "gpw", "normal": "npw", "tangent": "tpw", "interface": "ipw"}[family]] = 1.0
    _assert_parity(_data(family), weights)


@pytest.mark.parametrize("term", ALL_REGULARISATION)
def test_native_matches_cg_for_regularisation_term(term):
    # One term has weight 1. The other terms have the small background weight.
    kwargs = _background_regularisation(**{term: 1.0})
    _assert_parity(_data("value"), {**ZERO_DATA_WEIGHTS, "cpw": 1.0, **kwargs})


def test_native_matches_cg_for_border_regularisation():
    kwargs = _background_regularisation(dxx=1.0, dyy=1.0, dzz=1.0, dx=1.0, dy=1.0, dz=1.0)
    _assert_parity(_data("value"), {**ZERO_DATA_WEIGHTS, "cpw": 1.0, **kwargs})


def _vertical(points: np.ndarray) -> np.ndarray:
    return np.tile([0.0, 0.0, 1.0], (points.shape[0], 1))


def test_native_matches_cg_for_directional_regularisation():
    regularisation = {
        "isotropic": BACKGROUND,
        "directional": [{"weight": 1.0, "direction": _vertical, "name": "vertical smoothing"}],
    }
    _assert_parity(
        _data("value"),
        {**ZERO_DATA_WEIGHTS, "cpw": 1.0, "regularisation": regularisation},
    )


def test_native_node_order_matches_matrix_solver():
    # Field f = x + 2y + 3z. Value constraints at every node fix the field,
    # so interpolator.c must equal f at the nodes in the node order of the grid.
    n = 8
    grid = _make_grid(n)
    nodes = np.asarray(grid.nodes, dtype=float)
    expected = nodes[:, 0] + 2.0 * nodes[:, 1] + 3.0 * nodes[:, 2]
    data = np.column_stack([nodes, expected, np.full(len(nodes), 1.0)])
    setup = {**ZERO_DATA_WEIGHTS, "cpw": 1.0, **_background_regularisation()}
    # Weak regularisation: the data fixes the field.
    setup.update({name: 1e-4 for name in ALL_REGULARISATION})

    results = {}
    for solver in ("cg", "native"):
        builder = InterpolatorBuilder(
            interpolatortype=InterpolatorType.FINITE_DIFFERENCE,
            support=grid,
        )
        if solver == "native":
            builder.use_solver("native", **NATIVE_KWARGS)
        else:
            builder.use_solver(solver, **CG_KWARGS)
        builder.add_value_constraints(data).setup_interpolator(**setup).solve()
        results[solver] = np.asarray(builder.build().c, dtype=float).copy()

    span = float(np.ptp(expected))
    assert np.max(np.abs(results["cg"] - expected)) <= 1e-3 * span
    assert np.max(np.abs(results["native"] - results["cg"])) <= REL_TOL * span


def _bad_row_data() -> dict:
    """Good rows of every family plus rows that the matrix solver drops or tolerates."""
    data = _data("value")
    rng = np.random.default_rng(SEED + 1)
    value = data["add_value_constraints"]
    bad_value = np.array(
        [
            [10.0, 10.0, 10.0, 5.0, -1.0],  # negative weight
            [12.0, 11.0, 10.0, np.nan, 1.0],  # nan value
            [np.nan, 11.0, 10.0, 5.0, 1.0],  # nan position
        ]
    )
    p = _points(rng, 30)
    g = _field_gradient(p)
    g /= np.linalg.norm(g, axis=1)[:, None]
    bad_vectors = np.array([[0.0, 0.0, 0.0], [np.nan, 1.0, 0.0]])
    bad_p = _points(rng, 2)
    gradient = np.vstack([np.column_stack([p, g]), np.column_stack([bad_p, bad_vectors])])
    normal_p = _points(rng, 30)
    n = _field_gradient(normal_p)
    normal = np.column_stack([normal_p, n, np.ones(len(normal_p))])
    bad_normal = np.column_stack(
        [_points(rng, 3), np.array([[np.nan, 0, 1], [0, 0, 0], [0.1, 0.2, 0.3]]), [1.0, 1.0, -1.0]]
    )
    t = np.cross(_field_gradient(p), np.tile([0.0, 0.0, 1.0], (len(p), 1)))
    t /= np.linalg.norm(t, axis=1)[:, None]
    tangent = np.vstack([np.column_stack([p, t]), np.column_stack([bad_p, bad_vectors])])
    return {
        "add_value_constraints": np.vstack([value, bad_value]),
        "add_gradient_constraints": gradient,
        "add_normal_constraints": np.vstack([normal, bad_normal]),
        "add_tangent_constraints": tangent,
    }


def test_native_matches_cg_with_bad_data_rows():
    # Bad rows: negative weight, nan vector, nan value, nan position, zero vector.
    # The native solver must not raise. The field must match the matrix solver.
    weights = {**ZERO_DATA_WEIGHTS, "cpw": 1.0, "gpw": 1.0, "npw": 1.0, "tpw": 1.0}
    reference = _build("cg", _bad_row_data(), weights)
    native = _build("native", _bad_row_data(), weights)
    assert np.all(np.isfinite(reference))
    assert np.all(np.isfinite(native))
    span = float(np.ptp(reference))
    assert np.max(np.abs(native - reference)) <= REL_TOL * span


def _bad_interface_data() -> dict:
    """Interface groups of 3 points plus a nan id, a nan position and a one point group."""
    data = _data("interface")
    rng = np.random.default_rng(SEED + 2)
    extra = np.column_stack([_points(rng, 4), [np.nan, 0.0, 9.0, 4.0]])
    extra[1, 0] = np.nan  # nan position inside group 0
    # Group 4 keeps one valid point only: the other point has a nan position.
    extra[3, 1] = np.nan
    extra = np.vstack(
        [
            extra,
            np.column_stack([_points(rng, 1), [4.0]]),
            np.column_stack([_points(rng, 1), [9.0]]),
        ]
    )
    # The setter may drop some of these rows. Inject the rows after the setter.
    data["_interface_rows"] = np.vstack([data.pop("set_interface_constraints"), extra])
    return data


def _build_interface(solver: str, data: dict, ipw) -> np.ndarray:
    data = dict(data)
    rows = data.pop("_interface_rows")
    builder = InterpolatorBuilder(
        interpolatortype=InterpolatorType.FINITE_DIFFERENCE, support=_make_grid()
    )
    if solver == "native":
        builder.use_solver("native", **NATIVE_KWARGS)
    else:
        builder.use_solver(solver, **CG_KWARGS)
    for method, array in data.items():
        getattr(builder, method)(array)
    # Inject the rows so that they reach the adapter and the matrix path as they are.
    builder.interpolator.get_interface_constraints = lambda: rows
    builder.setup_interpolator(
        interface_pair_mode="star", **{**ZERO_DATA_WEIGHTS, "cpw": 1.0, "ipw": ipw}
    )
    builder.solve()
    return np.asarray(builder.build().c, dtype=float).copy()


@pytest.mark.parametrize("ipw", [1.0, -1.0])
def test_native_matches_cg_with_bad_interface_rows(ipw):
    # Bad rows: nan id, nan position, groups with one valid point.
    # A negative ipw acts like its absolute value on the matrix path.
    data = _bad_interface_data()
    reference = _build_interface("cg", data, ipw)
    native = _build_interface("native", data, ipw)
    assert np.all(np.isfinite(native))
    span = float(np.ptp(reference))
    assert np.max(np.abs(native - reference)) <= REL_TOL * span


def test_native_negative_interface_weight_equals_absolute_weight():
    data = _bad_interface_data()
    positive = _build_interface("native", data, 1.0)
    negative = _build_interface("native", data, -1.0)
    assert np.max(np.abs(positive - negative)) <= REL_TOL * float(np.ptp(positive))


def test_native_interface_only_single_point_groups_adds_nothing():
    # Every group has one point. Nothing reaches the core.
    from loop_interpolation import _native

    calls = []

    class _Problem:
        def add_interface(self, *args):
            calls.append(args)

    class _Interp:
        def __init__(self):
            self.interpolation_weights = {"ipw": 1.0}

        def __getattr__(self, name):
            return lambda: None

        def get_interface_constraints(self):
            return np.array([[1.0, 1, 1, 0], [2, 2, 2, 1], [3, 3, 3, np.nan], [4, 4, 4, np.nan]])

    _native._add_data(_Problem(), _Interp())
    assert calls == []


def test_native_mixed_interface_weights_still_raise():
    from loop_interpolation import _native

    class _Interp:
        def __init__(self):
            self.interpolation_weights = {"ipw": np.array([1.0, 2.0, 1.0])}

        def __getattr__(self, name):
            return lambda: None

        def get_interface_constraints(self):
            return np.array([[1.0, 1, 1, 0], [2, 2, 2, 0], [3, 3, 3, 0]])

    with pytest.raises(NotImplementedError):
        _native._add_data(object(), _Interp())
