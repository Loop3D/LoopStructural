"""Optional native and Surfe extension cells in a rotated frame (FR-030, SC-018, OP21).

These tests use the real extensions. They do not stub them. Each test skips when its
extension is not installed. The cells are ``native-extension`` and ``surfe-extension``.

The fit tolerances are not measured yet, because the extensions are absent in the
development environment (OP21). A limit is never relaxed. A miss is a fit-accuracy
finding.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.backend_matrix import require_optional
from tests.fixtures.coordinate_frames import (
    ANCHORS,
    QUERY_OFFSETS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    model_with_frame,
    r1_gradient,
    r1_interpolator,
    r1_observations,
    r1_value,
)

ANCHOR_IDS = [f"anchor{index}" for index in range(len(ANCHORS))]
ROTATION_IDS = [f"rot{int(angle)}" for angle in ROTATIONS_DEG]


def _assert_world_field(interp, anchor, observations) -> None:
    """Check the value and gradient of the R1 field at the two query points (SC-002)."""
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(observations["val"]).max())
    value_error = np.abs(interp.evaluate_value(points) - r1_value(points, anchor)).max()
    assert value_error <= 1e-8 * scale, f"world value error {value_error:.3e}"
    expected = r1_gradient(points, anchor)
    gradient = interp.evaluate_gradient(points)
    assert gradient.shape == expected.shape
    gradient_error = np.abs(gradient - expected)
    limit = 1e-8 * np.maximum(1.0, np.abs(expected))
    assert np.all(gradient_error <= limit), f"world gradient error {gradient_error.max():.3e}"


@pytest.mark.parametrize("rotation", ROTATIONS_DEG, ids=ROTATION_IDS)
@pytest.mark.parametrize("anchor", ANCHORS, ids=ANCHOR_IDS)
def test_native_extension_r1_fit_in_world(anchor, rotation):
    """The native solver fits the R1 field in a rotated frame and agrees with cg."""
    require_optional("native")
    # Fit support: the 27 value rows are projected to the local frame by the interpolator.
    interp, observations = r1_interpolator("fdi", anchor, rotation, exact=False)
    assert observations.shape[0] == 27
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(observations["val"]).max())
    # The helper solved the system with cg. Keep that field as the reference.
    cg_field = np.array(interp.c, copy=True)
    cg_values = interp.evaluate_value(points)
    assert np.all(np.isfinite(cg_field))

    assert interp.solve_system(solver="native", solver_kwargs={"tol": 1e-12})
    assert interp.up_to_date
    native_field = np.asarray(interp.c)
    assert native_field.shape == cg_field.shape
    assert np.all(np.isfinite(native_field))

    _assert_world_field(interp, anchor, observations)
    native_values = interp.evaluate_value(points)
    assert np.abs(native_values - cg_values).max() <= 1e-6 * scale
    assert np.abs(native_field - cg_field).max() <= 1e-6 * scale


@pytest.mark.parametrize("rotation", ROTATIONS_DEG, ids=ROTATION_IDS)
@pytest.mark.parametrize("anchor", ANCHORS, ids=ANCHOR_IDS)
def test_surfe_extension_r1_values_and_gradients_in_world(anchor, rotation):
    """The Surfe wrapper returns the R1 values and gradients in world coordinates."""
    require_optional("surfe")
    from loop_interpolation import SurfeRBFInterpolator

    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
    interp = SurfeRBFInterpolator()
    interp.bounding_box = model.bounding_box
    observations = r1_observations(anchor)
    rows = observations[["X", "Y", "Z", "val", "w"]].to_numpy()
    assert rows.shape[0] == 27
    interp.set_value_constraints(rows)
    # The default kernel stays. poly_order=1 lets the interpolant hold a linear trend.
    interp.setup_interpolator(poly_order=1)
    interp.solve_system()

    _assert_world_field(interp, anchor, observations)
    # The gradient constant must differ from zero, so that a lost rotation shows.
    assert np.abs(R1_GRADIENT).max() > 0.0
