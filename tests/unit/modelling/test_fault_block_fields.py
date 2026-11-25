"""Block fields: a feature evaluated as it is in one fault block.

The models are built in the unit cube. The stratigraphic data are made the
way LoopStructural evaluates a faulted feature: each random present-day point
is restored through the faults, and the analytic field (z) of the restored
point is its value. The interpolated field is then z in restored space, and
the faults are the only source of a jump.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import (
    OLD,
    YOUNG,
    add_fault,
    add_layers,
    dip_vectors,
    interior_points,
    plane_points,
    young_forward,
    young_side,
)


@pytest.fixture(scope="module")
def planar():
    """Planar normal fault, dip 60 toward +x, constant slip 0.15 down dip."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15)
    add_layers(model, [fault])
    return model


@pytest.fixture(scope="module")
def finite():
    """Finite fault (BaseFault3D) with a tip loop inside the model."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(
        model, "f1", (0.5, 0.5, 0.5), 70, 90, 0.2,
        profile="BaseFault3D", axes=(0.35, 0.3, 0.5), steps=10,
    )
    add_layers(model, [fault])
    return model


@pytest.fixture(scope="module")
def displaced():
    """An old fault displaced by a young fault (both constant slip)."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", YOUNG["center"], YOUNG["dip"], YOUNG["dip_azimuth"],
                       YOUNG["displacement"])
    plane, _, _ = plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=15)
    present = young_forward(plane)
    keep = np.abs(young_side(plane)) > 0.04
    center = young_forward(np.asarray([OLD["center"]]))[0]
    add_fault(model, "old", center, OLD["dip"], OLD["dip_azimuth"], OLD["displacement"],
               points=present[keep], faults=[young])
    return model


def test_planar_block_values_differ_by_the_analytic_offset(planar):
    pts = interior_points()
    hanging = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    footwall = planar.evaluate_feature_block_value("strat", pts, {"f1": -1})
    # Restoring the hanging wall moves a point 0.15 up dip: z grows by 0.15 sin(60).
    offset = 0.15 * np.sin(np.radians(60))
    np.testing.assert_allclose(hanging - footwall, offset, atol=2e-3)


def test_planar_block_value_equals_the_field_inside_its_block(planar):
    pts = interior_points()
    fault = planar["f1"]
    gx = fault.evaluate_frame_value(pts)
    value = planar.evaluate_feature_value("strat", pts)
    for side in (1, -1):
        inside = side * gx > 1e-6
        block = planar.evaluate_feature_block_value("strat", pts, {"f1": side})
        np.testing.assert_array_equal(block[inside], value[inside])


def test_planar_block_gradient_is_smooth_across_the_fault(planar):
    pts = interior_points(n=100)
    for side in (1, -1):
        grad = planar.evaluate_feature_block_gradient("strat", pts, {"f1": side})
        np.testing.assert_allclose(grad, np.tile([0.0, 0.0, 1.0], (len(pts), 1)), atol=2e-2)


def test_block_absent_fault_does_not_move_points(planar):
    pts = interior_points()
    absent = planar.evaluate_feature_block_value("strat", pts, {"f1": 0})
    np.testing.assert_allclose(absent, pts[:, 2], atol=2e-3)


def test_block_calls_leave_no_state(planar):
    pts = interior_points()
    before = planar.evaluate_feature_value("strat", pts)
    first = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    planar.evaluate_feature_block_value("strat", pts, {"f1": -1})
    again = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    np.testing.assert_array_equal(first, again)
    np.testing.assert_array_equal(before, planar.evaluate_feature_value("strat", pts))


def test_block_calls_in_threads_match_sequential_calls(planar):
    pts = interior_points(n=2000)
    blocks = [{"f1": 1}, {"f1": -1}] * 4
    expected = [planar.evaluate_feature_block_value("strat", pts, b) for b in blocks]
    with ThreadPoolExecutor(max_workers=4) as pool:
        got = list(pool.map(lambda b: planar.evaluate_feature_block_value("strat", pts, b), blocks))
    for e, g in zip(expected, got):
        np.testing.assert_array_equal(e, g)


def test_block_is_checked(planar):
    pts = interior_points(n=5)
    with pytest.raises(ValueError, match="not in the model"):
        planar.evaluate_feature_block_value("strat", pts, {"nofault": 1})
    with pytest.raises(ValueError, match="must be"):
        planar.evaluate_feature_block_value("strat", pts, {"f1": 2})
    with pytest.raises(TypeError):
        planar.evaluate_feature_block_value("strat", pts, [("f1", 1)])


def test_fault_object_can_be_a_block_key(planar):
    pts = interior_points(n=50)
    by_name = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    by_fault = planar.evaluate_feature_block_value("strat", pts, {planar["f1"]: 1})
    np.testing.assert_array_equal(by_name, by_fault)


def test_frame_value_is_a_signed_distance_near_the_fault(planar):
    fault = planar["f1"]
    normal, _, _ = dip_vectors(60, 90)
    pts = interior_points(n=200, lo=0.3, hi=0.7)
    distance = (pts - 0.5) @ normal
    gx = fault.evaluate_frame_value(pts)
    grad = fault.evaluate_frame_gradient(pts)
    signed = gx / np.linalg.norm(grad, axis=1)
    near = np.abs(distance) < 0.1
    np.testing.assert_allclose(signed[near], distance[near], atol=5e-3)


def test_constant_fault_has_no_tip(planar):
    pts = interior_points(n=10)
    assert np.all(np.isinf(planar["f1"].evaluate_tip_value(pts)))
    np.testing.assert_array_equal(planar["f1"].evaluate_tip_gradient(pts), 0.0)


def test_finite_fault_block_fields_are_equal_beyond_the_tip(finite):
    pts = np.random.default_rng(11).uniform(0.0, 1.0, size=(4000, 3))
    fault = finite["f1"]
    tip = fault.evaluate_tip_value(pts)
    hanging = finite.evaluate_feature_block_value("strat", pts, {"f1": 1})
    footwall = finite.evaluate_feature_block_value("strat", pts, {"f1": -1})
    beyond = tip < 0
    assert beyond.sum() > 100
    np.testing.assert_array_equal(hanging[beyond], footwall[beyond])
    # Inside the tip loop and near the fault the blocks differ.
    gx = fault.evaluate_frame_value(pts)
    active = (tip > 0.5) & (np.abs(gx) < 0.1)
    assert active.sum() > 10
    assert np.all(np.abs(hanging[active] - footwall[active]) > 1e-2)


def test_finite_fault_displacement_goes_to_zero_at_the_tip(finite):
    fault = finite["f1"]
    normal, down_dip, strike = dip_vectors(70, 90)
    # A line on the fault along strike, from the centre past the tip.
    s = np.linspace(0.0, 0.48, 241)
    line = np.asarray([0.5, 0.5, 0.5]) + s[:, None] * strike
    tip = fault.evaluate_tip_value(line)
    displacement = np.abs(fault.evaluate_displacement(line, block={"f1": 1}))
    assert tip[0] > 0.5 and tip[-1] < 0
    np.testing.assert_array_equal(displacement[tip < 0], 0.0)
    # Continuous: small just inside the tip.
    near = (tip > 0) & (tip < 0.02)
    assert near.any()
    assert np.all(displacement[near] < 0.01 * displacement.max())


def test_finite_fault_tip_gradient_matches_finite_differences(finite):
    fault = finite["f1"]
    pts = np.random.default_rng(5).uniform(0.3, 0.7, size=(50, 3))
    grad = fault.evaluate_tip_gradient(pts)
    h = 1e-4
    fd = np.zeros_like(pts)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        fd[:, axis] = (
            fault.evaluate_tip_value(pts + step) - fault.evaluate_tip_value(pts - step)
        ) / (2 * h)
    # Skip points near a kink of min(1 - |gy|, 1 - |gz|): |gy| = |gz|, gy = 0 or gz = 0.
    gy = np.abs(fault.evaluate_frame_value(pts, 1))
    gz = np.abs(fault.evaluate_frame_value(pts, 2))
    smooth = (np.abs(gy - gz) > 1e-2) & (gy > 1e-2) & (gz > 1e-2)
    assert smooth.sum() > 40
    np.testing.assert_allclose(grad[smooth], fd[smooth], atol=1e-6)


def test_displaced_fault_block_frames_match_the_undisplaced_fault(displaced):
    old = displaced["old"]
    plane, _, _ = plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=25, half=0.5)
    _, down_dip, _ = dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    # Footwall block of young: the old fault where it was. Hanging wall block:
    # the old fault moved forward by the slip of young, on both sides of young.
    for side, points in ((-1, plane), (1, plane + YOUNG["displacement"] * down_dip)):
        points = points[np.all((points > 0.05) & (points < 0.95), axis=1)]
        gx = old.evaluate_frame_value(points, block={"young": side})
        grad = old.evaluate_frame_gradient(points, block={"young": side})
        distance = gx / np.linalg.norm(grad, axis=1)
        assert np.nanmax(np.abs(distance)) < 1e-2


@pytest.fixture(scope="module")
def abutting():
    """A fault that abuts on the hanging wall of a main fault (a Y junction)."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    main = add_fault(model, "main", (0.4, 0.5, 0.5), 60, 90, 0.1)
    splay = add_fault(model, "splay", (0.7, 0.5, 0.6), 60, 270, 0.08)
    splay.add_abutting_fault(main, positive=True)
    add_layers(model, [splay, main])
    return model


def test_abutting_fault_is_absent_in_the_footwall_block_of_its_host(abutting):
    pts = interior_points(n=2000, lo=0.05, hi=0.95)
    values = {
        (m, s): abutting.evaluate_feature_block_value("strat", pts, {"main": m, "splay": s})
        for m in (1, -1)
        for s in (1, -1)
    }
    # In the footwall block of main the splay does not exist: its sides are equal.
    np.testing.assert_array_equal(values[(-1, 1)], values[(-1, -1)])
    # In the hanging wall block of main the splay offsets the layers.
    offset = np.abs(values[(1, 1)] - values[(1, -1)])
    assert np.median(offset) > 0.05
    # Absent is the same as the footwall of main for the splay.
    absent = abutting.evaluate_feature_block_value("strat", pts, {"main": -1, "splay": 0})
    np.testing.assert_array_equal(values[(-1, 1)], absent)


def test_displaced_fault_block_frame_equals_the_frame_on_its_side(displaced):
    old = displaced["old"]
    pts = np.random.default_rng(9).uniform(0.05, 0.95, size=(2000, 3))
    gx = old.evaluate_frame_value(pts)
    young_gx = displaced["young"].evaluate_frame_value(pts)
    for side in (1, -1):
        inside = side * young_gx > 1e-6
        block_gx = old.evaluate_frame_value(pts, block={"young": side})
        np.testing.assert_array_equal(block_gx[inside], gx[inside])


def test_restoring_points_on_a_younger_fault_gives_finite_points(displaced):
    # Points on the young fault get a NaN young-fault gradient; the old fault must not spread it.
    on_young, _, _ = plane_points(YOUNG["center"], YOUNG["dip"], YOUNG["dip_azimuth"], n=7, half=0.3)
    restored = displaced["old"].apply_to_points(on_young)
    assert np.isfinite(restored).all()
