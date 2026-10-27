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
CG_KWARGS = {"rtol": 1e-10, "atol": 1e-10, "maxiter": 50000}

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
        builder.use_solver("native")
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
            builder.use_solver("native")
        else:
            builder.use_solver(solver, **CG_KWARGS)
        builder.add_value_constraints(data).setup_interpolator(**setup).solve()
        results[solver] = np.asarray(builder.build().c, dtype=float).copy()

    span = float(np.ptp(expected))
    assert np.max(np.abs(results["cg"] - expected)) <= 1e-3 * span
    assert np.max(np.abs(results["native"] - results["cg"])) <= REL_TOL * span
