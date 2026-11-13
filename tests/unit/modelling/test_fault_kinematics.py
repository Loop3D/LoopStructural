"""Kinematics of fault restoration: forward movement and displaced faults."""

import numpy as np
import pytest

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import add_fault, dip_vectors, interior_points


@pytest.fixture(scope="module")
def planar_fault():
    """Planar normal fault, dip 60 toward +x, constant slip 0.15 down dip."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    return add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15)


def test_reverse_moves_the_hanging_wall_forward(planar_fault):
    _, down_dip, _ = dip_vectors(60, 90)
    pts = interior_points(n=200, lo=0.3, hi=0.7)
    hanging = planar_fault.evaluate_frame_value(pts) > 0.05
    moved = planar_fault.apply_to_points(pts, reverse=True)
    # Forward: the hanging wall moves 0.15 down dip; the footwall stays.
    expected = np.tile(0.15 * down_dip, (int(hanging.sum()), 1))
    np.testing.assert_allclose(moved[hanging] - pts[hanging], expected, atol=2e-3)
    footwall = planar_fault.evaluate_frame_value(pts) < 0
    np.testing.assert_array_equal(moved[footwall], pts[footwall])


def test_reverse_then_restore_returns_the_points(planar_fault):
    pts = interior_points(n=200, lo=0.3, hi=0.7)
    hanging = planar_fault.evaluate_frame_value(pts) > 0.05
    back = planar_fault.apply_to_points(planar_fault.apply_to_points(pts, reverse=True))
    np.testing.assert_allclose(back[hanging], pts[hanging], atol=2e-3)
