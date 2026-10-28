"""Tests for the errors of ``solver="native"`` (spec 001, task T029).

Two rules are checked here:

1. A setting that the native solver does not support raises
   ``NotImplementedError``. The message names the setting.
2. When the native core is not installed, choosing ``"native"`` raises
   ``ImportError`` with the install hint ``pip install loop-interpolation[native]``.

These tests do not need ``loop_interpolation_core``. The unsupported-setting
checks run before the core import. The missing-core check blocks the import
with ``sys.modules``.
"""

import re
import sys

import numpy as np
import pytest
from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid


def _make_grid() -> StructuredGrid:
    return StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([5, 5, 5]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )


def _region_mask(interp: FiniteDifferenceInterpolator) -> np.ndarray:
    # A boolean mask over all support nodes. All nodes are active.
    return np.ones(interp.support.n_nodes, dtype=bool)


def _build_interpolator() -> FiniteDifferenceInterpolator:
    interp = FiniteDifferenceInterpolator(_make_grid())
    # Five columns: X, Y, Z, value, weight.
    interp.set_value_constraints(
        np.array(
            [
                [1.0, 1.0, 1.0, 0.0, 1.0],
                [2.0, 2.0, 2.0, 1.0, 1.0],
                [3.0, 2.0, 1.0, 0.5, 1.0],
            ]
        )
    )
    return interp


# Each case gives:
#   setup_kwargs   - keywords for setup_interpolator. A callable value
#                    receives the interpolator and returns the real value.
#   attributes     - attributes set on the interpolator before setup.
#   solver_kwargs  - keywords for solve_system.
# Column scaling (apply_scaling_matrix) is not in this list. The native solver
# ignores it for any value (see the Answers in plan.md).
UNSUPPORTED_SETTINGS = [
    pytest.param(
        "region",
        {"setup_kwargs": {"region": _region_mask}},
        id="region",
    ),
    pytest.param(
        "use_regularisation_weight_scale",
        {"setup_kwargs": {"use_regularisation_weight_scale": True}},
        id="use_regularisation_weight_scale",
    ),
    pytest.param(
        "constant_norm_iterations",
        {"solver_kwargs": {"constant_norm_iterations": 3}},
        id="constant_norm_iterations",
    ),
    pytest.param(
        "interface_pair_mode",
        {"setup_kwargs": {"interface_pair_mode": "all"}},
        id="interface_pair_mode",
    ),
    pytest.param(
        "regularisation_matrix_free",
        {
            "setup_kwargs": {"regularisation_matrix_free": True},
            # The matrix-free flag is dropped when apply_scaling_matrix is True.
            "attributes": {"apply_scaling_matrix": False},
        },
        id="regularisation_matrix_free",
    ),
]


@pytest.mark.parametrize("setting, config", UNSUPPORTED_SETTINGS)
def test_unsupported_setting_with_native_raises_not_implemented(setting, config):
    interp = _build_interpolator()

    for name, value in config.get("attributes", {}).items():
        setattr(interp, name, value)

    setup_kwargs = {
        key: (value(interp) if callable(value) else value)
        for key, value in config.get("setup_kwargs", {}).items()
    }
    interp.setup_interpolator(**setup_kwargs)

    solver_kwargs = dict(config.get("solver_kwargs", {}))

    # The message must name the unsupported setting, so the user can find it.
    with pytest.raises(NotImplementedError, match=re.escape(setting)):
        interp.solve_system(solver="native", solver_kwargs=solver_kwargs)


def test_native_without_core_raises_import_error_with_install_hint(monkeypatch):
    # A None entry in sys.modules makes "import loop_interpolation_core" raise
    # ImportError. This works even when the core is installed.
    monkeypatch.setitem(sys.modules, "loop_interpolation_core", None)

    interp = _build_interpolator()
    interp.setup_interpolator()

    with pytest.raises(ImportError, match=re.escape("pip install loop-interpolation[native]")):
        interp.solve_system(solver="native")


# --- Interpolator classes and sparse-only rows (review finding 1) ---


def _fold_grid_data():
    return np.array([[1.0, 1.0, 1.0, 0.0, 1.0], [2.0, 2.0, 2.0, 1.0, 1.0]])


def test_native_rejects_p1_interpolator():
    from loop_common.supports import TetMesh
    from loop_interpolation import PiecewiseLinearInterpolator as P1

    interp = P1(TetMesh())
    with pytest.raises(NotImplementedError, match="FiniteDifferenceInterpolator"):
        interp.solve_system(solver="native")


def test_native_rejects_fd_fold_subclass():
    from loop_interpolation import FDFoldInterpolator

    interp = FDFoldInterpolator(_make_grid())
    with pytest.raises(NotImplementedError, match="FDFoldInterpolator"):
        interp.solve_system(solver="native")


def test_native_accepts_subclass_with_flag():
    from loop_interpolation import _native

    class Flagged(FiniteDifferenceInterpolator):
        supports_native = True

    _native.check_supported(Flagged(_make_grid()), {}, {})
    _native.check_supported(FiniteDifferenceInterpolator(_make_grid()), {}, {})


def test_native_rejects_subclass_without_flag():
    class Plain(FiniteDifferenceInterpolator):
        pass

    with pytest.raises(NotImplementedError, match="Plain"):
        Plain(_make_grid()).solve_system(solver="native")


def test_native_rejects_subclass_of_flagged_subclass():
    class Flagged(FiniteDifferenceInterpolator):
        supports_native = True

    class Child(Flagged):
        pass

    with pytest.raises(NotImplementedError, match="Child"):
        Child(_make_grid()).solve_system(solver="native")


def test_native_rejects_sparse_only_least_squares_rows():
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator()
    interp.add_constraints_to_least_squares(
        np.array([[1.0]]), np.array([0.0]), np.array([[0]]), name="extra"
    )
    with pytest.raises(NotImplementedError, match="add_constraints_to_least_squares"):
        interp.solve_system(solver="native")


def test_native_rejects_extra_rows_after_matrix_setup():
    interp = _build_interpolator()
    interp.solver = "cg"
    interp.setup_interpolator()
    assert len(interp.constraints) > 0
    interp.add_constraints_to_least_squares(
        np.array([[1.0]]), np.array([0.0]), np.array([[0]]), name="extra"
    )
    with pytest.raises(NotImplementedError, match="add_constraints_to_least_squares"):
        interp.solve_system(solver="native")


def test_native_setup_then_matrix_solver_raises_value_error():
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator()
    with pytest.raises(ValueError, match="setup"):
        interp.solve_system(solver="cg")
    # The failed call must not destroy the old field.
    assert not np.all(np.isnan(interp.c))


def test_native_setup_keeps_setup_kwargs_and_builds_no_rows():
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator(cpw=2.0)
    assert interp._setup_kwargs == {"cpw": 2.0}
    assert interp._native_setup
    assert len(interp.constraints) == 0
    assert len(interp.ineq_constraints) == 0
    assert interp._setup_constraint_names == frozenset()


def test_second_setup_gives_new_snapshot_and_flag():
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator()
    assert interp._native_setup
    # A new setup with a matrix solver clears the flag and records the new rows.
    interp.solver = "cg"
    interp.setup_interpolator()
    assert not interp._native_setup
    assert interp._setup_constraint_names == frozenset(interp.constraints)
    assert len(interp._setup_constraint_names) > 0
    assert interp.solve_system(solver="cg") in (True, False)
    # reset() clears the flag and the snapshot.
    interp.solver = "native"
    interp.setup_interpolator()
    interp.reset()
    assert not interp._native_setup
    assert interp._setup_constraint_names is None


def test_native_rejects_inequality_feature():
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator()
    interp.add_inequality_feature(lambda xyz: xyz[:, 2], lower=True)
    with pytest.raises(NotImplementedError, match="add_inequality_feature"):
        interp.solve_system(solver="native")


def test_native_accepts_normal_data_rows_after_native_setup():
    pytest.importorskip("loop_interpolation_core")
    interp = _build_interpolator()
    interp.solver = "native"
    interp.setup_interpolator()
    assert interp.solve_system(solver="native") in (True, False)
