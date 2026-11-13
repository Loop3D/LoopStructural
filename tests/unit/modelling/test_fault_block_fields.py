"""Block fields: a feature evaluated as it is in one fault block.

The models are built in the unit cube. The stratigraphic data are made the
way LoopStructural evaluates a faulted feature: each random present-day point
is restored through the faults, and the analytic field (z) of the restored
point is its value. The interpolated field is then z in restored space, and
the faults are the only source of a jump.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_function import (
    Composite,
    FaultDisplacement,
    Ones,
    Zeros,
)

COLUMNS = ["X", "Y", "Z", "val", "nx", "ny", "nz", "coord", "feature_name"]


def _dip_vectors(dip_deg, dip_azimuth_deg):
    """Unit normal (up), down-dip and strike vectors of a plane."""
    d, a = np.radians(dip_deg), np.radians(dip_azimuth_deg)
    horizontal = np.array([np.sin(a), np.cos(a), 0.0])
    down_dip = np.cos(d) * horizontal + np.array([0.0, 0.0, -np.sin(d)])
    normal = np.sin(d) * horizontal + np.array([0.0, 0.0, np.cos(d)])
    return normal, down_dip, np.cross(normal, down_dip)


def _constant_profile():
    """The hanging wall moves rigidly, the footwall does not move."""
    return FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)


def _plane_points(center, dip, dip_azimuth, n=11, half=0.9):
    normal, down_dip, strike = _dip_vectors(dip, dip_azimuth)
    u, v = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n), indexing="ij")
    pts = (
        np.asarray(center)
        + u.ravel()[:, None] * half * strike
        + v.ravel()[:, None] * half * down_dip
    )
    return pts, normal, down_dip


def _add_fault(model, name, center, dip, dip_azimuth, displacement, *, points=None,
               profile=None, axes=(4.0, 4.0, 4.0), faults=None, steps=1, nelements=4000):
    plane, normal, down_dip = _plane_points(center, dip, dip_azimuth)
    if points is None:
        points = plane
    inside = np.all((points > -0.1) & (points < 1.1), axis=1)
    points = points[inside]
    df = pd.DataFrame(
        {
            "X": points[:, 0], "Y": points[:, 1], "Z": points[:, 2], "val": 0.0,
            "nx": normal[0], "ny": normal[1], "nz": normal[2], "coord": 0,
            "feature_name": name,
        }
    )[COLUMNS]
    fault = model.create_and_add_fault(
        name,
        displacement,
        data=df,
        fault_slip_vector=-down_dip,  # restoration direction: normal fault
        fault_normal_vector=normal,
        fault_center=np.asarray(center, dtype=float),
        major_axis=axes[0],
        intermediate_axis=axes[1],
        minor_axis=axes[2],
        faultfunction=profile if profile is not None else _constant_profile(),
        faults=faults or [],
        nelements=nelements,
        fault_buffer=0.5,
    )
    fault.steps = steps
    return fault


def _add_layers(model, faults_oldest_first, n=800, seed=7, nelements=8000):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0.0, 1.0, size=(n, 3))
    restored = pts.copy()
    for f in reversed(faults_oldest_first):
        restored = f.apply_to_points(restored)
    df = pd.DataFrame(
        {
            "X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": restored[:, 2],
            "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": "strat",
        }
    )[COLUMNS]
    return model.create_and_add_foliation(
        "strat", data=df, faults=faults_oldest_first, nelements=nelements, buffer=0.3
    )


def _interior_points(n=400, seed=3, lo=0.15, hi=0.85):
    return np.random.default_rng(seed).uniform(lo, hi, size=(n, 3))


@pytest.fixture(scope="module")
def planar():
    """Planar normal fault, dip 60 toward +x, constant slip 0.15 down dip."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = _add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15)
    _add_layers(model, [fault])
    return model


@pytest.fixture(scope="module")
def finite():
    """Finite fault (BaseFault3D) with a tip loop inside the model."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = _add_fault(
        model, "f1", (0.5, 0.5, 0.5), 70, 90, 0.2,
        profile="BaseFault3D", axes=(0.35, 0.3, 0.5), steps=10,
    )
    _add_layers(model, [fault])
    return model


YOUNG = dict(center=(0.65, 0.5, 0.5), dip=65, dip_azimuth=270, displacement=0.12)
OLD = dict(center=(0.3, 0.5, 0.4), dip=55, dip_azimuth=90, displacement=0.1)


def _young_side(points):
    normal, _, _ = _dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    return (np.asarray(points) - np.asarray(YOUNG["center"])) @ normal


def _young_forward(points):
    """Move restored points forward through the young fault (exact: a translation)."""
    _, down_dip, _ = _dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    points = np.array(points, dtype=float)
    hanging = _young_side(points) > 0
    points[hanging] += YOUNG["displacement"] * down_dip
    return points


@pytest.fixture(scope="module")
def displaced():
    """An old fault displaced by a young fault (both constant slip)."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = _add_fault(model, "young", YOUNG["center"], YOUNG["dip"], YOUNG["dip_azimuth"],
                       YOUNG["displacement"])
    plane, _, _ = _plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=15)
    present = _young_forward(plane)
    keep = np.abs(_young_side(plane)) > 0.04
    center = _young_forward(np.asarray([OLD["center"]]))[0]
    _add_fault(model, "old", center, OLD["dip"], OLD["dip_azimuth"], OLD["displacement"],
               points=present[keep], faults=[young])
    return model


def test_planar_block_values_differ_by_the_analytic_offset(planar):
    pts = _interior_points()
    hanging = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    footwall = planar.evaluate_feature_block_value("strat", pts, {"f1": -1})
    # Restoring the hanging wall moves a point 0.15 up dip: z grows by 0.15 sin(60).
    offset = 0.15 * np.sin(np.radians(60))
    np.testing.assert_allclose(hanging - footwall, offset, atol=2e-3)


def test_planar_block_value_equals_the_field_inside_its_block(planar):
    pts = _interior_points()
    fault = planar["f1"]
    gx = fault.evaluate_frame_value(pts)
    value = planar.evaluate_feature_value("strat", pts)
    for side in (1, -1):
        inside = side * gx > 1e-6
        block = planar.evaluate_feature_block_value("strat", pts, {"f1": side})
        np.testing.assert_array_equal(block[inside], value[inside])


def test_planar_block_gradient_is_smooth_across_the_fault(planar):
    pts = _interior_points(n=100)
    for side in (1, -1):
        grad = planar.evaluate_feature_block_gradient("strat", pts, {"f1": side})
        np.testing.assert_allclose(grad, np.tile([0.0, 0.0, 1.0], (len(pts), 1)), atol=2e-2)


def test_block_absent_fault_does_not_move_points(planar):
    pts = _interior_points()
    absent = planar.evaluate_feature_block_value("strat", pts, {"f1": 0})
    np.testing.assert_allclose(absent, pts[:, 2], atol=2e-3)


def test_block_calls_leave_no_state(planar):
    pts = _interior_points()
    before = planar.evaluate_feature_value("strat", pts)
    first = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    planar.evaluate_feature_block_value("strat", pts, {"f1": -1})
    again = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    np.testing.assert_array_equal(first, again)
    np.testing.assert_array_equal(before, planar.evaluate_feature_value("strat", pts))


def test_block_calls_in_threads_match_sequential_calls(planar):
    pts = _interior_points(n=2000)
    blocks = [{"f1": 1}, {"f1": -1}] * 4
    expected = [planar.evaluate_feature_block_value("strat", pts, b) for b in blocks]
    with ThreadPoolExecutor(max_workers=4) as pool:
        got = list(pool.map(lambda b: planar.evaluate_feature_block_value("strat", pts, b), blocks))
    for e, g in zip(expected, got):
        np.testing.assert_array_equal(e, g)


def test_block_is_checked(planar):
    pts = _interior_points(n=5)
    with pytest.raises(ValueError, match="not in the model"):
        planar.evaluate_feature_block_value("strat", pts, {"nofault": 1})
    with pytest.raises(ValueError, match="must be"):
        planar.evaluate_feature_block_value("strat", pts, {"f1": 2})
    with pytest.raises(TypeError):
        planar.evaluate_feature_block_value("strat", pts, [("f1", 1)])


def test_fault_object_can_be_a_block_key(planar):
    pts = _interior_points(n=50)
    by_name = planar.evaluate_feature_block_value("strat", pts, {"f1": 1})
    by_fault = planar.evaluate_feature_block_value("strat", pts, {planar["f1"]: 1})
    np.testing.assert_array_equal(by_name, by_fault)


def test_frame_value_is_a_signed_distance_near_the_fault(planar):
    fault = planar["f1"]
    normal, _, _ = _dip_vectors(60, 90)
    pts = _interior_points(n=200, lo=0.3, hi=0.7)
    distance = (pts - 0.5) @ normal
    gx = fault.evaluate_frame_value(pts)
    grad = fault.evaluate_frame_gradient(pts)
    signed = gx / np.linalg.norm(grad, axis=1)
    near = np.abs(distance) < 0.1
    np.testing.assert_allclose(signed[near], distance[near], atol=5e-3)


def test_constant_fault_has_no_tip(planar):
    pts = _interior_points(n=10)
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
    normal, down_dip, strike = _dip_vectors(70, 90)
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
    plane, _, _ = _plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=25, half=0.5)
    _, down_dip, _ = _dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
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
    main = _add_fault(model, "main", (0.4, 0.5, 0.5), 60, 90, 0.1)
    splay = _add_fault(model, "splay", (0.7, 0.5, 0.6), 60, 270, 0.08)
    splay.add_abutting_fault(main, positive=True)
    _add_layers(model, [splay, main])
    return model


def test_abutting_fault_is_absent_in_the_footwall_block_of_its_host(abutting):
    pts = _interior_points(n=2000, lo=0.05, hi=0.95)
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
