"""The CG stopping tolerances that solve_with_cg fills in."""

import pytest

from loop_interpolation import _solver_strategy
from loop_interpolation._solver_strategy import DEFAULT_CG_RTOL, _set_cg_tolerances


@pytest.mark.parametrize(
    "given, tol, expected",
    [
        # A tol without an rtol: scipy's own rtol (1e-5) must not decide the stop.
        ({}, 1e-10, {"atol": 1e-10, "rtol": DEFAULT_CG_RTOL}),
        ({}, None, {"rtol": DEFAULT_CG_RTOL}),
        ({"rtol": 1e-3}, 1e-10, {"atol": 1e-10, "rtol": 1e-3}),
        ({"atol": 1e-6}, 1e-10, {"atol": 1e-6, "rtol": DEFAULT_CG_RTOL}),
        ({"atol": 1e-6}, None, {"atol": 1e-6}),
        ({"atol": 1e-6, "rtol": 1e-4}, 1e-10, {"atol": 1e-6, "rtol": 1e-4}),
    ],
)
def test_cg_tolerances(given, tol, expected):
    kwargs = dict(given)
    _set_cg_tolerances(kwargs, tol)
    assert kwargs == expected


def test_solve_with_cg_uses_them():
    assert _solver_strategy.solve_with_cg.__code__.co_names.count("_set_cg_tolerances") == 1
