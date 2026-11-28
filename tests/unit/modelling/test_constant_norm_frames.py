"""Rotated frame checks for the constant-norm FDI and PLI interpolators (R1, SC-002)."""

import itertools

import numpy as np
import pytest
from loop_common.supports import StructuredGrid, TetMesh
from loop_interpolation import ConstantNormFDIInterpolator, ConstantNormP1Interpolator

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    model_with_frame,
    r1_observations,
    r1_value,
)

NORM_TARGET = np.sqrt(0.29)
SOLVER_KWARGS = {"rtol": 1e-10, "atol": 1e-12, "maxiter": 20_000}
QUERY_OFFSETS = np.array([[2.2, 3.7, 4.3], [6.1, 8.2, 9.1]])

ENGINES = {
    "cn_fdi": (ConstantNormFDIInterpolator, StructuredGrid),
    "cn_pli": (ConstantNormP1Interpolator, TetMesh),
}


def _fit(engine, anchor, rotation):
    """Fit the R1 field with a constant-norm interpolator in a rotated frame."""
    interpolator_class, support_class = ENGINES[engine]
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
    bounding_box = model.bounding_box
    corners = np.array(list(itertools.product(*zip(bounding_box.origin, bounding_box.maximum))))
    local = bounding_box.project(corners)
    lower, upper = local.min(axis=0), local.max(axis=0)
    # 6 x 6 x 6 cells give the 216 requested cells.
    support = support_class(origin=lower, nsteps_cells=[6, 6, 6], step_vector=(upper - lower) / 6)
    interpolator = interpolator_class(support)
    interpolator.bounding_box = bounding_box
    interpolator.norm_length = NORM_TARGET
    observations = r1_observations(anchor)
    interpolator.set_value_constraints(observations[["X", "Y", "Z", "val", "w"]].to_numpy())
    interpolator.setup()
    assert interpolator.solve_system(solver="cg", solver_kwargs=SOLVER_KWARGS)
    return interpolator, observations


def _check(engine, anchor, rotation):
    interpolator, observations = _fit(engine, anchor, rotation)
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(observations["val"]).max())
    value_error = np.abs(interpolator.evaluate_value(points) - r1_value(points, anchor)).max()
    gradient_error = np.abs(interpolator.evaluate_gradient(points) - R1_GRADIENT).max()
    assert value_error <= 1e-8 * scale
    assert gradient_error <= 1e-8 * max(1.0, np.abs(R1_GRADIENT).max())


@pytest.mark.parametrize("engine", list(ENGINES))
def test_unrotated_anchor_zero_control(engine):
    """The unrotated frame at the origin is the control for the rotated cells."""
    _check(engine, ANCHORS[0], 0.0)


@pytest.mark.parametrize("engine", list(ENGINES))
@pytest.mark.parametrize("rotation", ROTATIONS_DEG)
@pytest.mark.parametrize("anchor_index", range(len(ANCHORS)))
def test_constant_norm_fit_recovers_r1_in_rotated_frame(engine, rotation, anchor_index):
    _check(engine, ANCHORS[anchor_index], rotation)
