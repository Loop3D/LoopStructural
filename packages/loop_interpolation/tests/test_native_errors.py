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
UNSUPPORTED_SETTINGS = [
    pytest.param(
        "region",
        {"setup_kwargs": {"region": _region_mask}},
        id="region",
    ),
    pytest.param(
        "apply_scaling_matrix",
        {"attributes": {"apply_scaling_matrix": True}},
        id="apply_scaling_matrix",
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
