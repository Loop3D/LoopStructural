"""Copy, stale frame and rebuild checks for the frame of an interpolator (U3.3, U3.9, SC-016)."""

import copy
import itertools

import numpy as np
import pytest
from loop_common.supports import P2TetMesh
from loop_interpolation import P2FoldInterpolator, StaleFrameError

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    LIFECYCLE_ENGINES,
    QUERY_OFFSETS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    AnalyticFold,
    model_with_frame,
    r1_interpolator,
    r1_observations,
    r1_value,
    r2_probes,
    rot_z,
    tp,
)
from tests.unit.modelling.test_fault_frame_analytic import _fault

COPY_ENGINES = ("fdi", "p2", "p2_fold")
POINTS_OFFSETS = np.array([[1.5, 2.5, 3.5], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])


def _p2_fold_interpolator(anchor, rotation_deg):
    """Return a P2 fold interpolator in a rotated frame, with no solve."""
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)
    bounding_box = model.bounding_box
    corners = np.array(list(itertools.product(*zip(bounding_box.origin, bounding_box.maximum))))
    local = bounding_box.project(corners)
    lower, upper = local.min(axis=0), local.max(axis=0)
    cells = (5, 6, 7)
    mesh = P2TetMesh(origin=lower, nsteps=cells, step_vector=(upper - lower) / cells)
    fold = AnalyticFold.constant(
        np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])
    )
    interp = P2FoldInterpolator(mesh, fold=fold)
    interp.bounding_box = bounding_box
    observations = r1_observations(anchor)
    interp.set_value_constraints(observations[["X", "Y", "Z", "val", "w"]].to_numpy())
    return interp


def _build(engine, anchor, rotation_deg):
    if engine == "p2_fold":
        return _p2_fold_interpolator(anchor, rotation_deg)
    return r1_interpolator(engine, anchor, rotation_deg, exact=True)[0]


def _stored_rows(interp):
    return interp.n_g + interp.n_i + interp.n_n + interp.n_t


@pytest.mark.parametrize("rotation", ROTATIONS_DEG, ids=lambda value: f"rot{value:g}")
@pytest.mark.parametrize("anchor_index", range(len(ANCHORS)), ids=lambda value: f"anchor{value}")
@pytest.mark.parametrize("engine", COPY_ENGINES)
def test_copy_keeps_frame(engine, anchor_index, rotation):
    """A copy shares the frame of its source and holds no constraint rows."""
    anchor = ANCHORS[anchor_index]
    interp = _build(engine, anchor, rotation)
    assert _stored_rows(interp) > 0
    duplicate = interp.copy()
    assert duplicate.bounding_box is interp.bounding_box
    points = anchor + POINTS_OFFSETS
    local = interp.world_to_local_points(points)
    vectors = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.3, -0.4, 0.5]])
    for method, data in (
        ("world_to_local_points", points),
        ("local_to_world_points", local),
        ("world_to_local_vectors", vectors),
        ("local_to_world_vectors", vectors),
    ):
        np.testing.assert_array_equal(
            getattr(duplicate, method)(data), getattr(interp, method)(data)
        )
    assert _stored_rows(duplicate) == 0
    assert duplicate.support is interp.support
    if engine == "p2_fold":
        assert duplicate.fold is interp.fold


@pytest.mark.parametrize("engine", LIFECYCLE_ENGINES)
def test_stale_frame_rejected_without_state_change(engine):
    """A changed transform makes each use raise once and leaves the state intact."""
    anchor = ANCHORS[1]
    interp, _ = r1_interpolator(engine, anchor, 0.0, exact=False)
    bounding_box = interp.bounding_box
    points = anchor + QUERY_OFFSETS
    snapshot = {
        "data": copy.deepcopy(interp.data),
        "c": copy.deepcopy(interp.c),
        "up_to_date": interp.up_to_date,
        "valid": interp.valid,
        "counts": (interp.n_g, interp.n_i, interp.n_n, interp.n_t),
        "report": interp.latest_diagnostics_report,
    }
    bounding_box.set_local_transform(local_origin=anchor + (3, -2, 1), rotation_matrix=rot_z(30.0))
    rows = r1_observations(anchor)[["X", "Y", "Z", "val", "w"]].to_numpy()
    assert rows.shape[0] == 27
    calls = (
        lambda: interp.evaluate_value(points),
        lambda: interp.evaluate_gradient(points),
        lambda: interp.set_value_constraints(rows),
        lambda: interp.to_dict(),
    )
    for call in calls:
        with pytest.raises(StaleFrameError, match="^the local transform changed"):
            call()
    assert interp.data.keys() == snapshot["data"].keys()
    for key, value in snapshot["data"].items():
        np.testing.assert_array_equal(interp.data[key], value)
    np.testing.assert_array_equal(interp.c, snapshot["c"])
    assert interp.up_to_date == snapshot["up_to_date"]
    assert interp.valid == snapshot["valid"]
    assert (interp.n_g, interp.n_i, interp.n_n, interp.n_t) == snapshot["counts"]
    assert interp.latest_diagnostics_report is snapshot["report"]


def test_restored_transform_is_usable():
    """Setting the old transform back makes the interpolator usable again."""
    anchor = ANCHORS[1]
    interp, _ = r1_interpolator("fdi", anchor, 0.0, exact=False)
    bounding_box = interp.bounding_box
    points = anchor + QUERY_OFFSETS
    first = interp.evaluate_value(points)
    bounding_box.set_local_transform(local_origin=anchor + (3, -2, 1), rotation_matrix=rot_z(30.0))
    bounding_box.set_local_transform(local_origin=anchor, rotation_matrix=rot_z(0.0))
    np.testing.assert_array_equal(interp.evaluate_value(points), first)


@pytest.mark.parametrize("engine", ("fdi", "pli", "p2", "cn_fdi", "cn_pli"))
def test_rebuild_in_new_frame_preserves_world_evaluation(engine):
    """A new interpolator in the new frame gives the same world field (SC-002)."""
    anchor = ANCHORS[1]
    # P2 uses the exact nodal control: its fit misses SC-002 in every frame (OP7).
    interp, observations = r1_interpolator(engine, anchor, 30.0, exact=(engine == "p2"))
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(observations["val"]).max())
    value_error = np.abs(interp.evaluate_value(points) - r1_value(points, anchor)).max()
    gradient_error = np.abs(interp.evaluate_gradient(points) - R1_GRADIENT).max()
    assert value_error <= 1e-8 * scale
    assert gradient_error <= 1e-8 * max(1.0, np.abs(R1_GRADIENT).max())


@pytest.mark.parametrize("engine", ("FDI", "PLI"))
def test_fault_restoration_uses_current_frame(engine):
    """A stale frame is refused when points are restored through a fault frame.

    No restoration cache exists in the code. The check is that a stale
    frame is refused and that a rebuild equals an uncached control.
    """
    anchor = ANCHORS[1]
    points = r2_probes(anchor)[:200]
    fault = _fault(anchor, 0.0, engine)
    # Restore once, so that a cache (if it existed) would hold this result.
    fault.apply_to_points(points)
    fault.builder.model.bounding_box.set_local_transform(
        local_origin=anchor + (3, -2, 1), rotation_matrix=rot_z(30.0)
    )
    with pytest.raises(StaleFrameError, match="^the local transform changed"):
        fault.apply_to_points(points)
    # Build the same fault in a fresh model with the new frame, and in a second
    # fresh model as the uncached control.
    rebuilt = _fault(anchor, 30.0, engine).apply_to_points(points)
    control = _fault(anchor, 30.0, engine).apply_to_points(points)
    tolerance = tp(20.0, np.abs(points).max())
    assert np.all(np.isfinite(rebuilt))
    assert np.max(np.abs(rebuilt - control)) <= tolerance
