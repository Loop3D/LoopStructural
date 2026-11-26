"""Gradient constraints and the per-point weight column, for each discrete engine."""

import logging

import numpy as np
import pytest
from loop_common.supports import P2TetMesh
from loop_interpolation import (
    FiniteDifferenceInterpolator,
    P1Interpolator,
    P2Interpolator,
    StructuredGrid,
    TetMesh,
)

GRID = {"origin": np.zeros(3), "nsteps": np.array([5, 5, 5]), "step_vector": np.ones(3)}
ENGINES = {
    "FDI": lambda: FiniteDifferenceInterpolator(StructuredGrid(**GRID)),
    "P1": lambda: P1Interpolator(TetMesh(**GRID)),
    "P2": lambda: P2Interpolator(P2TetMesh(**GRID)),
}
CENTRE = np.array([[2.5, 2.5, 2.5]])
# two values on a vertical line fix the sign; alone they give a vertical gradient
VALUES = np.array([[2.5, 2.5, 1.0, 0.0, 1.0], [2.5, 2.5, 4.0, 3.0, 1.0]])


def _gradient_at_centre(interpolator):
    interpolator.setup()
    interpolator.solve_system("cg")
    g = interpolator.support.evaluate_gradient(CENTRE, interpolator.c)[0]
    return g / np.linalg.norm(g)


def _gradient_row(vector, w):
    return np.hstack([CENTRE[0], vector, w])


@pytest.mark.parametrize("engine", ENGINES)
def test_gradient_constraints_set_the_direction(engine):
    interpolator = ENGINES[engine]()
    interpolator.set_value_constraints(VALUES)
    interpolator.set_gradient_constraints(_gradient_row([1.0, 0.0, 1.0], 1.0)[None, :])
    g = _gradient_at_centre(interpolator)
    np.testing.assert_allclose(g, [2**-0.5, 0.0, 2**-0.5], atol=0.1)


@pytest.mark.parametrize("engine", ENGINES)
def test_the_weight_column_is_used(engine):
    # value: two points at one place, the heavy one wins
    interpolator = ENGINES[engine]()
    interpolator.set_value_constraints(
        np.vstack([VALUES, [2.5, 2.5, 2.5, 0.0, 1.0], [2.5, 2.5, 2.5, 3.0, 100.0]])
    )
    interpolator.setup()
    interpolator.solve_system("cg")
    value = interpolator.support.evaluate_value(CENTRE, interpolator.c)[0]
    assert value > 2.5
    # gradient: two opposite tilts at one place, the heavy one wins
    interpolator = ENGINES[engine]()
    interpolator.set_value_constraints(VALUES)
    interpolator.set_gradient_constraints(
        np.vstack([_gradient_row([1.0, 0.0, 1.0], 1.0), _gradient_row([-1.0, 0.0, 1.0], 100.0)])
    )
    g = _gradient_at_centre(interpolator)
    assert g[0] < -0.5


@pytest.mark.parametrize(
    "engine, family, data",
    [
        ("P1", "interface", [[1.0, 1.0, 1.0, 0.0, 1.0], [3.0, 3.0, 3.0, 0.0, 1.0]]),
        ("P2", "interface", [[1.0, 1.0, 1.0, 0.0, 1.0], [3.0, 3.0, 3.0, 0.0, 1.0]]),
        ("P2", "inequality", [[1.0, 1.0, 1.0, 0.0, 1.0, 1.0]]),
        ("P2", "inequality_pairs", [[1.0, 1.0, 1.0, 0.0, 1.0], [3.0, 3.0, 3.0, 1.0, 1.0]]),
    ],
)
def test_ignored_data_gives_a_warning(engine, family, data, caplog):
    interpolator = ENGINES[engine]()
    interpolator.set_value_constraints(VALUES)
    interpolator.data[family] = np.array(data)
    with caplog.at_level(logging.WARNING):
        interpolator.setup()
    assert f"{family} constraints" in caplog.text
