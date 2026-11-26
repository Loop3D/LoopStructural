"""apply_to_points of a fault whose frame is cut by another fault.

On a fault that cuts the frame of another fault, the frame has no single
gradient, so evaluate_gradient gives NaN there (#421). The points must still
move, with the direction of the frame without the cutting faults. A frame
without cutting faults keeps the NaN (spec 005, test_fault_frame_validation).
"""

import numpy as np

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import OLD, YOUNG, add_fault, interior_points


def _nan_on_every_other_row(member):
    original = member.evaluate_gradient

    def gradient(points, *args, **kwargs):
        g = np.array(original(points, *args, **kwargs), dtype=float)
        g[::2] = np.nan
        return g

    return gradient


def _cut_pair():
    """The old fault is cut by the young fault, so its frame members have a fault."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", **YOUNG, steps=3)
    old = add_fault(model, "old", **OLD, steps=3, faults=[young])
    model.update()  # the frame members get the cutting fault when they are built
    assert old[1].faults
    return old


def test_nan_gradient_of_a_cut_frame_still_moves_the_points(monkeypatch):
    old = _cut_pair()
    points = interior_points()
    expected = old.apply_to_points(points.copy())
    moved = ~np.isclose(expected, points).all(axis=1)
    assert moved.any(), "the fault must move some points"

    monkeypatch.setattr(old[1], "evaluate_gradient", _nan_on_every_other_row(old[1]))
    restored = old.apply_to_points(points.copy())

    assert not np.isnan(restored).any()
    # The fallback direction is close to the cut-frame direction, so the points move
    # about as far as before.
    np.testing.assert_allclose(
        np.linalg.norm(restored - points, axis=1),
        np.linalg.norm(expected - points, axis=1),
        rtol=0.05,
        atol=1e-9,
    )


def test_point_with_a_nan_coordinate_is_not_moved(monkeypatch):
    old = _cut_pair()
    points = np.vstack((interior_points(n=20), [[np.nan, 0.5, 0.5]]))
    monkeypatch.setattr(old[1], "evaluate_gradient", _nan_on_every_other_row(old[1]))
    restored = old.apply_to_points(points.copy())
    np.testing.assert_array_equal(restored[-1], points[-1])
    assert not np.isnan(restored[:-1]).any()


def test_frame_without_cutting_faults_keeps_a_nan_gradient(monkeypatch):
    """Without cutting faults, a NaN frame gradient still gives NaN (no fallback)."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(model, "young", **YOUNG, steps=3)
    points = interior_points()
    moved = ~np.isclose(fault.apply_to_points(points.copy()), points).all(axis=1)
    monkeypatch.setattr(fault[1], "evaluate_gradient", _nan_on_every_other_row(fault[1]))
    restored = fault.apply_to_points(points.copy())
    assert np.isnan(restored[moved]).any()
