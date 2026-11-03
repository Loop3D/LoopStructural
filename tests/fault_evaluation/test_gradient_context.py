"""Tests for the sharing of fault restorations in gradient evaluation (D5, OP1, T020).

OP1 rule, as accepted for this spec: inside one top-level evaluation call, one
restoration is made for each distinct ORDERED fault set and each distinct point
batch. The recursive batches of fault frames are counted separately. The total
across different batches is NOT claimed to be F.

A gradient of a faulted feature uses two point batches:

* the N original points, for the region masks (the regions apply the faults of
  their own feature);
* the 4N corner points of the tetrahedra, for the faults of the feature itself.

The persistent cache of the model is OFF in all tests of this file. Only the
sharing inside one evaluation scope is measured.

Each count test records a diagnostics dictionary (calls, distinct batches, batch
sizes, restored points and restored bytes). The dictionary is printed and given
to ``record_property``. Tests that check numbers (golden gradients, mask pattern,
group priority, output array) pass before and after the sharing change. They
guard the numerical behaviour. Tests that check call counts fail while equal
batches are restored more than once.
"""

from __future__ import annotations

import collections
import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from LoopStructural.modelling.features._base_geological_feature import BaseFeature
from LoopStructural.modelling.features._evaluation import evaluation_scope, restore
from LoopStructural.modelling.features._geological_feature import GeologicalFeature
from LoopStructural.modelling.features.fault import FaultSegment
from LoopStructural.utils.maths import gradient_from_tetrahedron, regular_tetraherdron_for_points
from tests.fault_evaluation.golden import assert_matches
from tests.fault_evaluation.models import build_model

DATA = Path(__file__).parent / "data"
MANIFEST = json.loads((DATA / "manifest.json").read_text(encoding="utf8"))

# Counts that the accepted OP1 rule gives for the golden points (see the report of T020).
F4_U3_GROUP_GRADIENT_CALLS = 8  # 4 faults on N points + 4 faults on 4N points
F4_U3_GROUP_GRADIENT_FAULT_POINTS_FACTOR = 20  # 4 * N + 4 * 4N
CUTTING_FAULTS_GROUP_0_GRADIENT_CALLS = 25  # distinct (fault, batch) pairs, frames included


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _golden(name: str) -> dict[str, np.ndarray]:
    """Return the stored legacy reference of a golden model."""
    with np.load(DATA / f"{name}.npz", allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def _assert_golden(name: str, actual: dict[str, np.ndarray]) -> None:
    """Compare chosen outputs with the stored legacy output, with the golden tolerances."""
    expected = _golden(name)
    reference = {key: expected[key] for key in actual}
    assert_matches(name, dict(actual), reference, MANIFEST["models"][name]["box_size"])


def _model(name: str):
    """Build a golden model with the persistent cache off (only scope sharing is measured)."""
    model = build_model(name)
    model.evaluation_cache.enabled = False
    return model


def _group(model, index: int) -> GeologicalFeature:
    return model.get_feature_by_name(f"group_{index}")


def _digest(points: np.ndarray) -> bytes:
    return hashlib.md5(np.ascontiguousarray(points).tobytes()).digest()


class _Recorder:
    """Record every call of ``FaultSegment.apply_to_points``: fault, batch and depth."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, int, bytes, bool, int]] = []
        self._depth = 0
        self._original = FaultSegment.apply_to_points

    def __enter__(self) -> _Recorder:  # noqa: PYI034
        recorder, original = self, self._original

        def spy(fault, points, *args, **kwargs):
            reverse = bool(kwargs.get("reverse", args[0] if args else False))
            array = np.asarray(points)
            recorder.calls.append(
                (fault.name, len(array), _digest(array), reverse, recorder._depth)
            )
            recorder._depth += 1
            try:
                return original(fault, points, *args, **kwargs)
            finally:
                recorder._depth -= 1

        FaultSegment.apply_to_points = spy
        return self

    def __exit__(self, *exc_info) -> None:
        FaultSegment.apply_to_points = self._original

    @property
    def total(self) -> int:
        return len(self.calls)

    @property
    def duplicates(self) -> int:
        """Number of calls that repeat an earlier (fault, batch, reverse) call."""
        keys = [(name, digest, reverse) for name, _, digest, reverse, _ in self.calls]
        return len(keys) - len(set(keys))

    def batch_calls(self, digest: bytes) -> int:
        """Number of calls whose input is the batch with this digest."""
        return sum(1 for call in self.calls if call[2] == digest)

    def diagnostics(self) -> dict:
        """Return the batch, count and memory diagnostics of the recorded calls."""
        rows = sum(call[1] for call in self.calls)
        return {
            "fault_calls_total": self.total,
            "fault_calls_top_level": sum(1 for call in self.calls if call[4] == 0),
            "duplicate_calls": self.duplicates,
            "distinct_batches": len({call[2] for call in self.calls}),
            "calls_by_batch_size": dict(collections.Counter(call[1] for call in self.calls)),
            "restored_points": rows,
            # Memory per actual restored point: 3 float64 coordinates, corners included.
            "restored_bytes": rows * 3 * np.dtype(float).itemsize,
        }


def _record(label, recorder, record_property) -> dict:
    """Print the diagnostics of a recorder and give them to the test report."""
    diagnostics = recorder.diagnostics()
    print(f"[gradient diagnostics] {label}: {json.dumps(diagnostics, sort_keys=True)}")
    record_property(f"gradient_diagnostics_{label}", json.dumps(diagnostics, sort_keys=True))
    return diagnostics


def _default_scale(model) -> float:
    return float(np.min(model.bounding_box.step_vector) / 10)


def _corners(model, points: np.ndarray, scale: float | None = None) -> np.ndarray:
    """Return the 4N corner points that ``evaluate_gradient`` builds for ``points``."""
    scale = _default_scale(model) if scale is None else scale
    return np.array(regular_tetraherdron_for_points(points, scale).reshape(-1, 3))


def _final_scale(feature, points: np.ndarray, scale: float) -> float:
    """Replay the refinement loop of ``evaluate_gradient`` and return the final scale."""
    resolved = False
    tetrahedron = regular_tetraherdron_for_points(points, scale)
    while not resolved:
        resolved = True
        for fault in feature.faults:
            values = fault[0].evaluate_value(tetrahedron.reshape(-1, 3), fillnan="nearest")
            values = values.reshape(tetrahedron.shape[0], 4)
            if np.any(~(np.all(values > 0, axis=1) | np.all(values < 0, axis=1))):
                scale *= 0.5
                tetrahedron = regular_tetraherdron_for_points(points, scale)
                resolved = False
    return scale


def _final_corners(model, feature, points: np.ndarray) -> np.ndarray:
    """Return the 4N corners that the refinement loop gives for the default scale."""
    return _corners(model, points, _final_scale(feature, points, _default_scale(model)))


def _legacy_gradient(feature, points: np.ndarray, scale: float) -> np.ndarray:
    """Recompute a faulted gradient with the legacy steps and no scope.

    ``scale`` is the start scale. The refinement loop is replayed first.
    The mask uses the original points. The corners are restored by the ordered faults.
    """
    result = np.full(points.shape, np.nan)
    mask = feature._calculate_mask(points)
    tetrahedron = regular_tetraherdron_for_points(points, _final_scale(feature, points, scale))
    restored = np.array(tetrahedron.reshape(-1, 3))
    for fault in feature.faults if feature.faults_enabled else []:
        restored = fault.apply_to_points(restored)
    values = feature.interpolator.evaluate_value(restored).reshape(-1, 4)
    result[mask, :] = gradient_from_tetrahedron(tetrahedron[mask, :, :], values[mask])
    return result


def _legacy_model_gradient(model, points: np.ndarray) -> np.ndarray:
    """Combine the group gradients with the priority of ``evaluate_model_gradient``.

    Every group call is its own top-level call. The youngest valid group wins.
    """
    result = np.zeros(points.shape)
    for group in reversed(model.stratigraphic_column.get_groups()):
        feature = model.get_feature_by_name(group.name)
        gradient = feature.evaluate_gradient(np.array(points))
        valid = ~np.isnan(gradient).any(axis=1)
        result[valid] = gradient[valid]
    return result


# ---------------------------------------------------------------------------
# Counts with one fault set: f4_u3 (4 faults, 3 unconformities)
# ---------------------------------------------------------------------------


def test_group_gradient_restores_original_and_corner_batches_once(record_property):
    """``group_3`` gradient: one restoration per fault for N points and for 4N corners.

    The three unconformity regions use the same ordered set of 4 faults on the N
    original points as the mask of the group. Expected: 4 calls on N points and 4 calls
    on 4N points = 8. Legacy: 16 (the regions repeat the N batch three times).
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    n = len(points)
    with _Recorder() as recorder:
        gradient = _group(model, 3).evaluate_gradient(points)
    diagnostics = _record("f4_u3_group_3", recorder, record_property)

    _assert_golden("f4_u3", {"points": points, "feature::group_3::gradient": gradient})
    assert diagnostics["calls_by_batch_size"] == {n: 4, 4 * n: 4}
    assert diagnostics["duplicate_calls"] == 0
    assert diagnostics["fault_calls_total"] == F4_U3_GROUP_GRADIENT_CALLS
    assert diagnostics["restored_points"] == F4_U3_GROUP_GRADIENT_FAULT_POINTS_FACTOR * n
    assert diagnostics["restored_bytes"] == F4_U3_GROUP_GRADIENT_FAULT_POINTS_FACTOR * n * 24


def test_region_masks_and_corners_use_the_expected_batches():
    """The first fault of the ordered set gets the original points and the exact corners.

    The masks are made from the original points. They never use restored points. The corner
    batch is the 4N points of the tetrahedra of the original points.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    with _Recorder() as recorder:
        gradient = _group(model, 3).evaluate_gradient(points)
    original = _digest(points)
    corners = _digest(_final_corners(model, _group(model, 3), points))

    # Only the first fault of the set sees an unrestored batch.
    assert recorder.batch_calls(original) == 1
    assert recorder.batch_calls(corners) == 1
    first = _group(model, 3).faults[0].name
    assert [call[0] for call in recorder.calls if call[2] in (original, corners)] == [first, first]
    np.testing.assert_array_equal(
        np.isnan(gradient).any(axis=1), ~_group(model, 3)._calculate_mask(points)
    )


def test_model_gradient_restores_each_batch_once_for_all_groups(record_property):
    """``evaluate_model_gradient`` shares N and 4N batches over all 4 groups.

    The 4 groups and their regions use one ordered set of 4 faults. The N batch and
    the 4N batch are each restored once: 8 calls. Legacy: 52.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    n = len(points)
    with _Recorder() as recorder:
        gradient = model.evaluate_model_gradient(points)
    diagnostics = _record("f4_u3_model_gradient", recorder, record_property)

    assert diagnostics["calls_by_batch_size"] == {n: 4, 4 * n: 4}
    assert diagnostics["duplicate_calls"] == 0
    assert diagnostics["fault_calls_total"] == 8
    assert diagnostics["restored_points"] == 20 * n

    expected = _legacy_model_gradient(_model("f4_u3"), points)
    np.testing.assert_allclose(gradient, expected, rtol=0.0, atol=1e-12, equal_nan=True)


def test_model_gradient_group_priority_is_unchanged():
    """The youngest group with a valid gradient wins. This is a guard that passes today.

    The golden group gradients give the reference. The golden groups do not overlap, so the
    next test checks the priority with overlapping groups.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    gradient = model.evaluate_model_gradient(points)

    expected = np.zeros(points.shape)
    for group in reversed(model.stratigraphic_column.get_groups()):
        group_gradient = golden[f"feature::{group.name}::gradient"]
        valid = ~np.isnan(group_gradient).any(axis=1)
        expected[valid] = group_gradient[valid]
    np.testing.assert_allclose(gradient, expected, rtol=0.0, atol=1e-8 * 1000.0)


def test_model_gradient_priority_with_overlapping_groups():
    """Where two groups are valid, the youngest group (Group_0) gives the gradient.

    The younger group is valid for x > 0.5. The older group is valid everywhere and has
    another gradient. This guard passes before and after the sharing change.
    """
    from LoopStructural import GeologicalModel
    from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn

    class Field:
        def __init__(self, gradient):
            self.gradient = np.array(gradient)

        def evaluate_value(self, points):
            return points @ self.gradient

        def evaluate_gradient(self, points):
            return np.broadcast_to(self.gradient, points.shape).copy()

    model = GeologicalModel([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    model.evaluation_cache.enabled = False
    young = GeologicalFeature(
        "young", None, regions=[lambda p: p[:, 0] > 0.5], faults=[], interpolator=Field([1, 0, 0])
    )
    old = GeologicalFeature("old", None, faults=[], interpolator=Field([0, 0, 2]))
    model._add_feature(young)
    model._add_feature(old)
    column = StratigraphicColumn()
    column.clear(basement=False)
    for index in (1, 0):
        high = column.add_unit(f"unit_{index}_high")
        low = column.add_unit(f"unit_{index}_low")
        for unit, (low_value, high_value) in ((high, (5.0, 10.0)), (low, (0.0, 5.0))):
            unit.min_value, unit.max_value = low_value, high_value
        if index:
            column.add_unconformity("unconformity_1")
    column.group_mapping.update({"Group_0": "young", "Group_1": "old"})
    model.stratigraphic_column = column
    points = np.array([[0.25, 0.5, 0.5], [0.75, 0.5, 0.5]])

    gradient = model.evaluate_model_gradient(points)

    np.testing.assert_allclose(gradient, [[0.0, 0.0, 2.0], [1.0, 0.0, 0.0]], atol=1e-12)


def test_repeated_groups_in_one_scope_share_the_batches(record_property):
    """Two groups with the same ordered set, inside one scope: N and 4N batches once each."""
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    with _Recorder() as recorder, evaluation_scope(model):
        third = _group(model, 3).evaluate_gradient(points)
        second = _group(model, 2).evaluate_gradient(points)
        again = _group(model, 3).evaluate_gradient(points)
    diagnostics = _record("f4_u3_repeated_groups", recorder, record_property)

    assert diagnostics["fault_calls_total"] == 8
    assert diagnostics["duplicate_calls"] == 0
    _assert_golden(
        "f4_u3",
        {
            "points": points,
            "feature::group_3::gradient": third,
            "feature::group_2::gradient": second,
        },
    )
    np.testing.assert_array_equal(third, again)


# ---------------------------------------------------------------------------
# Element scale and refinement
# ---------------------------------------------------------------------------


def test_different_element_scales_make_different_corner_batches(record_property):
    """Equal scales share the 4N batch. A different scale is another batch.

    Three gradients in one scope: scales 20, 20 and 7. The N batch is restored once
    (4 calls). The 4N batch of scale 20 is restored once and the one of scale 7 once
    (8 calls). Total 12. Never merge two scales: the gradients must differ.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    n = len(points)
    with _Recorder() as recorder, evaluation_scope(model):
        first = feature.evaluate_gradient(points, element_scale_parameter=20.0)
        repeated = feature.evaluate_gradient(points, element_scale_parameter=20.0)
        other = feature.evaluate_gradient(points, element_scale_parameter=7.0)
    diagnostics = _record("f4_u3_scales", recorder, record_property)

    assert diagnostics["calls_by_batch_size"] == {n: 4, 4 * n: 8}
    assert diagnostics["duplicate_calls"] == 0
    assert diagnostics["fault_calls_total"] == 12

    reference = _model("f4_u3")
    np.testing.assert_array_equal(first, repeated)
    np.testing.assert_allclose(
        first,
        _legacy_gradient(_group(reference, 3), points, 20.0),
        rtol=0.0,
        atol=1e-10,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        other,
        _legacy_gradient(_group(reference, 3), points, 7.0),
        rtol=0.0,
        atol=1e-10,
        equal_nan=True,
    )
    assert not np.allclose(first, other, equal_nan=True, rtol=0.0, atol=1e-12)


def test_tetrahedron_refinement_restores_only_the_final_corners(caplog, record_property):
    """A refined tetrahedron is restored once, at the final scale. No stale batch is restored.

    At the default scale the corners of the golden points touch the fault frames, so the loop
    halves the scale (the log shows the warning). The final scale comes from an independent
    replay of the loop. Only the corners of the final scale are restored, and the gradient equals
    the legacy gradient at that scale.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    n = len(points)
    feature = _group(model, 3)
    start = _default_scale(model)
    final = _final_scale(feature, points, start)
    assert final < start, "these points must need a refinement"
    with caplog.at_level(logging.WARNING), _Recorder() as recorder:
        gradient = feature.evaluate_gradient(points)
    diagnostics = _record("f4_u3_refinement", recorder, record_property)

    assert any("Refining the tetrahedron" in message for message in caplog.messages)
    assert recorder.batch_calls(_digest(_corners(model, points, start))) == 0
    assert recorder.batch_calls(_digest(_corners(model, points, final))) == 1
    assert diagnostics["calls_by_batch_size"] == {n: 4, 4 * n: 4}
    assert diagnostics["fault_calls_total"] == 8
    np.testing.assert_allclose(
        gradient,
        _legacy_gradient(_group(_model("f4_u3"), 3), points, start),
        rtol=0.0,
        atol=1e-10,
        equal_nan=True,
    )


# ---------------------------------------------------------------------------
# Different fault sets
# ---------------------------------------------------------------------------


def test_reordered_fault_list_is_a_different_set_in_model_gradient(record_property):
    """group_1 has the faults in reverse order: set B. The other features use set A.

    A and B are two keys, for the N batch and for the 4N batch: (4 + 4) + (4 + 4) = 16.
    Legacy: 52. A merged key would give 8 and wrong values.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    group = _group(model, 1)
    group.faults = list(reversed(group.faults))
    assert group.faults != _group(model, 3).faults
    points = golden["points"]
    n = len(points)
    with _Recorder() as recorder:
        gradient = model.evaluate_model_gradient(points)
    diagnostics = _record("f4_u3_reordered_model_gradient", recorder, record_property)

    assert diagnostics["calls_by_batch_size"] == {n: 8, 4 * n: 8}
    assert diagnostics["duplicate_calls"] == 0
    assert diagnostics["fault_calls_total"] == 16

    reference = _model("f4_u3")
    reference_group = _group(reference, 1)
    reference_group.faults = list(reversed(reference_group.faults))
    np.testing.assert_allclose(
        gradient, _legacy_model_gradient(reference, points), rtol=0.0, atol=1e-12, equal_nan=True
    )
    # The order changes the restored corners, so a merged key would give wrong values.
    corners = _corners(model, points)
    ordered = reverse_ordered = corners
    for fault in _group(reference, 3).faults:
        ordered = fault.apply_to_points(ordered)
    for fault in reference_group.faults:
        reverse_ordered = fault.apply_to_points(reverse_ordered)
    assert not np.allclose(ordered, reverse_ordered, rtol=0.0, atol=1e-9)


def test_onlap_domain_and_stopped_fault_model_gradients(record_property):
    """Models with an empty fault set and with a shared set: counts and golden gradients.

    ``stopped_fault`` has the sets [] and [fault_0]: only the 4N batch of the group with
    the fault is restored (1 call). ``onlap_domain`` shares [fault_0] over 3 groups and 2
    regions: one call on N points and one on 4N points (2 calls). Legacy: 4.
    """
    for name, expected_calls in (("stopped_fault", 1), ("onlap_domain", 2)):
        model, golden = _model(name), _golden(name)
        points = golden["points"]
        with _Recorder() as recorder:
            gradient = model.evaluate_model_gradient(points)
        diagnostics = _record(f"{name}_model_gradient", recorder, record_property)
        assert diagnostics["fault_calls_total"] == expected_calls, name
        assert diagnostics["duplicate_calls"] == 0, name
        expected = _legacy_model_gradient(_model(name), points)
        np.testing.assert_allclose(gradient, expected, rtol=0.0, atol=1e-12, equal_nan=True)


# ---------------------------------------------------------------------------
# Child-fault frames
# ---------------------------------------------------------------------------


def test_child_fault_frame_batches_are_counted_separately(record_property):
    """``cutting_faults`` group_0 gradient: no duplicate (fault, batch) pair remains.

    group_0 uses [fault_1, fault_0]. The frame of fault_1 uses [fault_0] and restores its
    own nested batches, for N points and for the corners. Every nested batch counts
    on its own (OP1). With one restoration per distinct batch there are 25 calls
    (legacy: 46), and none repeats an earlier fault and batch.
    """
    model, golden = _model("cutting_faults"), _golden("cutting_faults")
    points = golden["points"]
    with _Recorder() as recorder:
        gradient = _group(model, 0).evaluate_gradient(points)
    diagnostics = _record("cutting_faults_group_0", recorder, record_property)

    _assert_golden("cutting_faults", {"points": points, "feature::group_0::gradient": gradient})
    assert diagnostics["duplicate_calls"] == 0
    assert diagnostics["fault_calls_total"] == CUTTING_FAULTS_GROUP_0_GRADIENT_CALLS
    # Nested calls come from the frame of fault_1 and are not top-level calls.
    assert diagnostics["fault_calls_top_level"] < diagnostics["fault_calls_total"]


# ---------------------------------------------------------------------------
# Reverse restoration, disabled faults
# ---------------------------------------------------------------------------


def test_reverse_and_forward_restoration_of_corners_are_two_batches(record_property):
    """Reverse and forward restoration of the same corners are two keys, in one scope.

    1 forward call and 1 reverse call = 2 (one fault). The results equal the manual restorations
    (the legacy reverse path gives the values of the forward path, but the keys stay separate).
    """
    model = _model("f1_u0")
    feature = _group(model, 0)
    # Two points near the fault, so the restoration moves the corners.
    near = np.array([[450.0, 500.0, -500.0], [550.0, 500.0, -450.0]])
    corners = _corners(model, near)
    with _Recorder() as recorder, evaluation_scope(model):
        forward = restore(feature, corners)
        backward = restore(feature, corners, reverse=True)
        forward_again = restore(feature, corners)
        backward_again = restore(feature, corners, reverse=True)
    _record("f1_u0_reverse_corners", recorder, record_property)

    expected_forward, expected_backward = corners, corners
    for fault in feature.faults:
        expected_forward = fault.apply_to_points(expected_forward)
    for fault in feature.faults:
        expected_backward = fault.apply_to_points(expected_backward, reverse=True)
    np.testing.assert_array_equal(forward, expected_forward)
    np.testing.assert_array_equal(backward, expected_backward)
    np.testing.assert_array_equal(forward, forward_again)
    np.testing.assert_array_equal(backward, backward_again)
    assert recorder.total == 2


def test_disabled_group_faults_leave_only_the_region_batches(record_property):
    """Only ``group_3`` has ``faults_enabled = False``. Its regions keep their faults.

    The group restores no corner batch. The regions restore the N batch once: 4 calls.
    The gradient equals the legacy gradient without faults. Legacy: 12 calls.
    """
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    feature.faults_enabled = False
    with _Recorder() as recorder:
        gradient = feature.evaluate_gradient(points)
    diagnostics = _record("f4_u3_group_disabled", recorder, record_property)

    assert diagnostics["calls_by_batch_size"] == {len(points): 4}
    assert diagnostics["fault_calls_total"] == 4
    expected = _legacy_gradient_no_faults(_group(_model("f4_u3"), 3), points)
    np.testing.assert_allclose(gradient, expected, rtol=0.0, atol=1e-10, equal_nan=True)


def _legacy_gradient_no_faults(feature, points: np.ndarray) -> np.ndarray:
    """Legacy gradient of a feature whose own faults are off (regions keep their faults)."""
    feature.faults_enabled = False
    scale = _default_scale(feature.model)
    return _legacy_gradient(feature, points, scale)


def test_all_faults_disabled_make_no_call():
    """With every fault switch off, no restoration is made. The result is NaN-compatible."""
    model, golden = _model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    seen: set[int] = set()

    def disable(item) -> None:
        if id(item) in seen or not isinstance(item, BaseFeature):
            return
        seen.add(id(item))
        item.faults_enabled = False
        for region in item.regions:
            disable(region)

    disable(feature)
    with _Recorder() as recorder:
        gradient = feature.evaluate_gradient(points)
    assert recorder.total == 0
    assert gradient.shape == points.shape
    np.testing.assert_array_equal(np.isnan(gradient).any(axis=1), ~feature._calculate_mask(points))


# ---------------------------------------------------------------------------
# Faulted fold frames
# ---------------------------------------------------------------------------


def test_faulted_fold_frame_gradients_share_the_corner_batch(record_property):
    """Fold frame coordinates, the folded unit and group_0 use [fault_0] on the same corners.

    Inside one scope all five gradients restore the 4N batch once: 1 call. Without
    sharing the five gradients restore it 6 times (the cross product evaluates two
    coordinates). The gradients equal the golden gradients.
    """
    model, golden = _model("faulted_fold"), _golden("faulted_fold")
    points = golden["points"]
    frame = model.get_feature_by_name("fold_frame")
    with _Recorder() as recorder, evaluation_scope(model):
        results = {
            "fold_frame__0": frame[0].evaluate_gradient(points),
            "fold_frame__1": frame[1].evaluate_gradient(points),
            "fold_frame__2": frame[2].evaluate_gradient(points),
            "folded_unit": model.get_feature_by_name("folded_unit").evaluate_gradient(points),
            "group_0": _group(model, 0).evaluate_gradient(points),
        }
    diagnostics = _record("faulted_fold_scope", recorder, record_property)

    _assert_golden(
        "faulted_fold",
        {"points": points}
        | {f"feature::{name}::gradient": gradient for name, gradient in results.items()},
    )
    assert diagnostics["calls_by_batch_size"] == {4 * len(points): 1}
    assert diagnostics["fault_calls_total"] == 1
    assert diagnostics["restored_points"] == 4 * len(points)


# ---------------------------------------------------------------------------
# Numerical guards that pass before and after the change
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["f4_u3", "cutting_faults", "faulted_fold", "onlap_domain"])
def test_model_gradient_matches_golden_group_priority(name):
    """The model gradient equals the golden group gradients with the youngest-valid priority."""
    model, golden = _model(name), _golden(name)
    points = golden["points"]
    gradient = model.evaluate_model_gradient(points)

    expected = np.zeros(points.shape)
    for group in reversed(model.stratigraphic_column.get_groups()):
        group_gradient = golden[f"feature::{group.name}::gradient"]
        valid = ~np.isnan(group_gradient).any(axis=1)
        expected[valid] = group_gradient[valid]
    box = MANIFEST["models"][name]["box_size"]
    np.testing.assert_allclose(gradient, expected, rtol=0.0, atol=1e-8 * box)


class _IdentityFault(BaseFeature):
    """A fault double that leaves coordinates unchanged."""

    def __init__(self, name="identity"):
        super().__init__(name)

    def evaluate_value(self, points, **kwargs):
        return points[:, 0]

    def evaluate_gradient(self, points, **kwargs):
        return np.broadcast_to([1.0, 0.0, 0.0], points.shape).copy()

    def __getitem__(self, key):
        assert key == 0
        return self

    def apply_to_points(self, points, reverse=False):
        return points

    def copy(self, name=None):
        return _IdentityFault(name or self.name)

    def get_data(self, value_map=None):
        return []


class _LinearField:
    """A linear interpolator with a constant gradient."""

    gradient = np.array([2.0, -3.0, 0.5])

    def evaluate_value(self, points):
        return points @ self.gradient + 7.0

    def evaluate_gradient(self, points):
        return np.broadcast_to(self.gradient, points.shape).copy()


def test_faulted_gradient_output_array_shape_and_fill():
    """D0 regression: the output has shape (N, 3) and a region-masked point stays NaN.

    This is an equivalent of ``tests/unit/modelling/test_faulted_gradient_output.py``.
    The sharing of batches must not change the shape, the dtype or the NaN fill.
    """
    points = np.array([[1.0, 2.0, 3.0], [-1.0, 2.0, 3.0], [4.0, -2.0, 1.0]])
    feature = GeologicalFeature(
        "linear",
        builder=None,
        regions=[lambda p: p[:, 0] > 0],
        faults=[_IdentityFault("a"), _IdentityFault("b")],
        interpolator=_LinearField(),
    )
    gradient = feature.evaluate_gradient(points, element_scale_parameter=1)

    assert gradient.shape == (3, 3)
    assert gradient.dtype == np.float64
    np.testing.assert_allclose(gradient[[0, 2]], np.tile(_LinearField.gradient, (2, 1)), atol=1e-12)
    assert np.isnan(gradient[1]).all()
