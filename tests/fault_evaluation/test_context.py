"""Tests for the evaluation scope and the sharing of fault restorations (D3, D5, OP1).

Two kinds of test are in this file.

* Behaviour tests use the public evaluation functions. They count the calls of
  ``FaultSegment.apply_to_points`` and compare the values with the stored golden
  output (or with an independent recomputation that does not use a scope).
  The golden models are small, so each test uses the stored golden points.
* Private API tests check ``LoopStructural.modelling.features._evaluation``.
  They fail with a clear message while this module is absent.

OP1 rule, as accepted for this spec: inside one top-level evaluation call, one
restoration is made for each distinct ORDERED fault set and each distinct point
batch. Recursive fault-frame batches are counted separately. The total count
across different batches is NOT claimed to be F.
"""

from __future__ import annotations

import importlib
import inspect
import json
from pathlib import Path

import numpy as np
import pytest

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation
from LoopStructural import GeologicalModel
from LoopStructural.modelling.features._base_geological_feature import BaseFeature
from LoopStructural.modelling.features._geological_feature import GeologicalFeature
from LoopStructural.modelling.features.fault import FaultSegment
from tests.fault_evaluation.golden import assert_matches
from tests.fault_evaluation.models import build_model

MODULE_NAME = "LoopStructural.modelling.features._evaluation"
DATA = Path(__file__).parent / "data"
MANIFEST = json.loads((DATA / "manifest.json").read_text(encoding="utf8"))

# The present signatures. The performance change must not alter them.
EXPECTED_SIGNATURES = {
    (GeologicalModel, "evaluate_model"): (
        "(self, xyz: 'np.ndarray', *, scale: 'bool' = True) -> 'np.ndarray'"
    ),
    (GeologicalModel, "evaluate_model_gradient"): (
        "(self, points: 'np.ndarray', *, scale: 'bool' = True) -> 'np.ndarray'"
    ),
    (GeologicalFeature, "evaluate_value"): (
        "(self, pos: 'np.ndarray', ignore_regions=False, fillnan=None) -> 'np.ndarray'"
    ),
    (GeologicalFeature, "evaluate_gradient"): (
        "(self, pos: 'np.ndarray', ignore_regions=False, element_scale_parameter=None)"
        " -> 'np.ndarray'"
    ),
    (BaseFeature, "evaluate_value"): "(self, pos, ignore_regions=False)",
    (BaseFeature, "evaluate_gradient"): "(self, pos, ignore_regions=False)",
}


def _api(name: str):
    """Return an object of the evaluation module, or fail with a clear reason when it is absent."""
    try:
        module = importlib.import_module(MODULE_NAME)
    except ImportError as error:
        pytest.fail(f"{MODULE_NAME} is required by D3 (import failed: {error})")
    value = getattr(module, name, None)
    assert value is not None, f"{MODULE_NAME}.{name} is required by D3"
    return value


def _golden(name: str) -> dict[str, np.ndarray]:
    """Return the stored legacy reference of a golden model."""
    with np.load(DATA / f"{name}.npz", allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def _assert_golden(name: str, actual: dict[str, np.ndarray]) -> None:
    """Compare chosen outputs with the stored legacy output, with the golden tolerances."""
    expected = _golden(name)
    actual = dict(actual)
    if "model_ids" in actual:
        # The comparator accepts a changed id only at a unit boundary.
        actual["unit_boundary_distances"] = expected["unit_boundary_distances"]
    reference = {key: expected[key] for key in actual}
    assert_matches(name, actual, reference, MANIFEST["models"][name]["box_size"])


def _counted(call):
    """Run ``call`` and return its result and a snapshot of the fault call counters.

    ``fault_calls_total`` is the nested total: every call of ``apply_to_points``,
    also the calls made inside the frame of another fault (OP1: counted separately).
    ``fault_calls_top_level`` counts only the calls that are not inside another
    ``apply_to_points`` call.
    """
    with EvaluationInstrumentation() as instrumentation:
        result = call()
    return result, instrumentation.snapshot()


def _group(model: GeologicalModel, index: int) -> GeologicalFeature:
    return model.get_feature_by_name(f"group_{index}")


def _reference_value(feature, points: np.ndarray, ignore_regions: bool = False) -> np.ndarray:
    """Recompute ``evaluate_value`` with the legacy steps and no scope.

    The mask uses the original points. The faults are applied one by one, in the
    order of ``feature.faults``.
    """
    mask = np.ones(len(points), dtype=bool)
    if not ignore_regions:
        for region in feature.regions:
            mask = np.logical_and(mask, region(points))
    restored = points
    if feature.faults_enabled:
        for fault in feature.faults:
            restored = fault.apply_to_points(restored)
    values = np.full(len(points), np.nan)
    values[mask] = feature.interpolator.evaluate_value(restored[mask, :])
    return values


def _manual_restore(feature, points: np.ndarray, reverse: bool = False) -> np.ndarray:
    for fault in feature.faults:
        points = fault.apply_to_points(points, reverse=reverse)
    return points


def _disable_faults(feature, seen: set[int] | None = None) -> None:
    """Switch off the faults of a feature and of all its region features."""
    seen = set() if seen is None else seen
    if id(feature) in seen or not isinstance(feature, BaseFeature):
        return
    seen.add(id(feature))
    feature.faults_enabled = False
    for region in feature.regions:
        _disable_faults(region, seen)


def _assert_no_stale_scope(feature, points: np.ndarray, model=None) -> None:
    """Show that no scope is active: a new scope must apply every fault of ``feature`` again."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    _, counts = _counted(lambda: _enter_and_restore(scope, restore, model, feature, points))
    assert counts["fault_calls_total"] == len(feature.faults), "a scope is still active"


def _enter_and_restore(scope, restore, model, feature, points):
    with scope(model):
        return restore(feature, points)


# ---------------------------------------------------------------------------
# Public signatures (these freeze the present behaviour and pass today)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("owner,name", list(EXPECTED_SIGNATURES))
def test_public_signature_is_unchanged(owner, name):
    """The public evaluation signatures are the same as before the change."""
    actual = str(inspect.signature(getattr(owner, name)))
    assert actual == EXPECTED_SIGNATURES[(owner, name)]


# ---------------------------------------------------------------------------
# One-call restoration sharing, model f4_u3 (4 faults, 3 unconformities)
# ---------------------------------------------------------------------------


def test_value_of_oldest_group_restores_once_per_fault():
    """The oldest group and its 3 unconformity regions use one ordered set of 4 faults.

    Legacy: 4 applications for the feature and 4 for each of 3 regions = 16.
    All applications are top level, because no fault frame is faulted.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    values, counts = _counted(lambda: _group(model, 3).evaluate_value(points))
    _assert_golden("f4_u3", {"points": points, "feature::group_3::value": values})
    assert counts["fault_calls_total"] == 4
    assert counts["fault_calls_top_level"] == 4


def test_evaluate_model_restores_once_per_fault():
    """``evaluate_model`` is one top-level call: 4 faults, once, for all 4 groups.

    Legacy: 52 applications (13 uses of the same ordered set of 4 faults).
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    ids, counts = _counted(lambda: model.evaluate_model(golden["points"]))
    _assert_golden("f4_u3", {"points": golden["points"], "model_ids": ids})
    assert counts["fault_calls_total"] == 4
    assert counts["fault_calls_top_level"] == 4


def test_gradient_values_stay_frozen():
    """The gradient algorithm is not changed by this task (T015 wraps value entry points only).

    The gradient of ``group_3`` and the model gradient equal the stored golden output.
    No count is checked here: the sharing of gradient batches is task T020/T021.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    gradient = _group(model, 3).evaluate_gradient(points)
    _assert_golden("f4_u3", {"points": points, "feature::group_3::gradient": gradient})
    model_gradient = model.evaluate_model_gradient(points)
    # The model gradient is the gradient of the youngest group that is not NaN.
    expected = np.zeros(points.shape)
    for group in reversed(model.stratigraphic_column.get_groups()):
        group_gradient = golden[f"feature::{group.name}::gradient"]
        valid = ~np.isnan(group_gradient).any(axis=1)
        expected[valid] = group_gradient[valid]
    np.testing.assert_allclose(model_gradient, expected, rtol=0.0, atol=1e-8 * 1000.0)


# ---------------------------------------------------------------------------
# Different and reordered fault sets
# ---------------------------------------------------------------------------


def test_stopped_fault_applies_each_distinct_set_once():
    """``stopped_fault`` has the sets [] (group_0, unconformity) and [fault_0] (group_1).

    The empty set needs no restoration. The set [fault_0] needs one. The sum over
    the distinct sets of the number of faults is 1 (legacy: also 1).
    """
    model, golden = build_model("stopped_fault"), _golden("stopped_fault")
    ids, counts = _counted(lambda: model.evaluate_model(golden["points"]))
    _assert_golden("stopped_fault", {"points": golden["points"], "model_ids": ids})
    assert counts["fault_calls_total"] == 1
    assert counts["fault_calls_top_level"] == 1


def test_onlap_domain_restores_the_shared_set_once():
    """In ``onlap_domain`` 3 groups and 2 regions use the ordered set [fault_0].

    Expected: 1 application. Legacy: 4.
    """
    model, golden = build_model("onlap_domain"), _golden("onlap_domain")
    ids, counts = _counted(lambda: model.evaluate_model(golden["points"]))
    _assert_golden("onlap_domain", {"points": golden["points"], "model_ids": ids})
    assert counts["fault_calls_total"] == 1


def test_reordered_fault_list_is_a_different_set():
    """Give group_1 the faults in the reverse order of the other features.

    group_1 then uses set B = [0, 1, 2, 3]. Its region __group_0_unconformity and
    all other features use set A = [3, 2, 1, 0]. Restoring in another order gives
    another result, so the two sets are two keys: 4 + 4 = 8 applications.
    An implementation that merges the sets would give 4 and wrong values.
    Legacy: 12 (group_1 and its region from group_1 use B, the other region uses A).
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    group = _group(model, 1)
    group.faults = list(reversed(group.faults))
    assert group.faults != _group(model, 3).faults
    values, counts = _counted(lambda: group.evaluate_value(points))
    np.testing.assert_allclose(values, _reference_value(group, points), atol=1e-8, equal_nan=True)
    assert counts["fault_calls_total"] == 8
    # The order changes the result, so the test would detect a merged key.
    assert not np.allclose(
        _manual_restore(group, points), _manual_restore(_group(model, 3), points)
    )


def _reordered_model() -> GeologicalModel:
    """Return ``f4_u3`` where group_1 has the faults in the reverse order (sets A and B)."""
    model = build_model("f4_u3")
    group = _group(model, 1)
    group.faults = list(reversed(group.faults))
    return model


def test_reordered_fault_list_in_evaluate_model():
    """The model call with sets A and B applies each ordered set once: 4 + 4 = 8."""
    model, golden = _reordered_model(), _golden("f4_u3")
    ids, counts = _counted(lambda: model.evaluate_model(golden["points"]))
    assert counts["fault_calls_total"] == 8
    assert ids.shape == (len(golden["points"]),)


# ---------------------------------------------------------------------------
# reverse, disabled faults, ignore_regions, zero displacement (values only)
# ---------------------------------------------------------------------------


def test_reverse_restoration_equals_the_legacy_result():
    """``_apply_faults(reverse=True)`` gives the legacy sequence of reversed restorations."""
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    actual, counts = _counted(lambda: feature._apply_faults(points, reverse=True))
    np.testing.assert_array_equal(actual, _manual_restore(feature, points, reverse=True))
    assert counts["fault_calls_total"] == 4
    assert not np.allclose(actual, points)


def test_disabled_faults_make_no_application():
    """With ``faults_enabled = False`` no fault is applied and the values do not change."""
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    _disable_faults(feature)
    values, counts = _counted(lambda: feature.evaluate_value(points))
    assert counts["fault_calls_total"] == 0
    np.testing.assert_allclose(values, _reference_value(feature, points), atol=1e-8, equal_nan=True)
    assert not np.allclose(values, golden["feature::group_3::value"], equal_nan=True)


def test_group_disabled_but_regions_enabled():
    """Only the group has ``faults_enabled = False``. Its regions keep the same ordered faults.

    The key includes ``faults_enabled``, so the group must not use the restoration of
    its regions (and the regions must not use the empty result of the group).
    Expected: 4 applications (the regions; the group applies none). Legacy: 12.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    feature.faults_enabled = False
    values, counts = _counted(lambda: feature.evaluate_value(points))
    np.testing.assert_allclose(values, _reference_value(feature, points), atol=1e-8, equal_nan=True)
    assert counts["fault_calls_total"] == 4


def test_restore_enabled_and_disabled_feature_in_one_scope():
    """Same points and same ordered faults, one feature enabled and one disabled, one scope.

    The enabled feature applies the faults (4). The disabled feature gets equal points and
    makes no new application. The result must not depend on the order of the two calls.
    """
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, disabled = golden["points"], _group(model, 3)
    enabled = next(
        region
        for region in disabled.regions
        if isinstance(region, BaseFeature) and region.faults == disabled.faults
    )
    expected = _manual_restore(enabled, points)
    disabled.faults_enabled = False
    try:

        def run(first, second):
            with scope(model):
                return restore(first, points), restore(second, points)

        (restored_on, restored_off), counts = _counted(lambda: run(enabled, disabled))
        assert counts["fault_calls_total"] == 4
        np.testing.assert_array_equal(restored_on, expected)
        np.testing.assert_array_equal(restored_off, points)
        (restored_off, restored_on), counts = _counted(lambda: run(disabled, enabled))
        assert counts["fault_calls_total"] == 4
        np.testing.assert_array_equal(restored_on, expected)
        np.testing.assert_array_equal(restored_off, points)
    finally:
        disabled.faults_enabled = True


def test_disabled_faults_in_evaluate_model():
    """The model call with all faults off applies no fault."""
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    for feature in model.features:
        _disable_faults(feature)
    ids, counts = _counted(lambda: model.evaluate_model(golden["points"]))
    assert counts["fault_calls_total"] == 0
    assert ids.shape == (len(golden["points"]),)


def test_ignore_regions_equals_the_legacy_result():
    """With ``ignore_regions=True`` the regions are not used: only the feature restores."""
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    values, counts = _counted(lambda: feature.evaluate_value(points, ignore_regions=True))
    np.testing.assert_allclose(
        values, _reference_value(feature, points, ignore_regions=True), atol=1e-8, equal_nan=True
    )
    assert counts["fault_calls_total"] == 4
    assert np.isfinite(values).all()


def test_zero_displacement_equals_the_unfaulted_result():
    """A fault with displacement 0 moves no point. The values equal the values without faults."""
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    for fault in _group(model, 3).faults:
        fault.displacement = 0.0
    values = _group(model, 3).evaluate_value(points)
    unfaulted = build_model("f4_u3")
    feature = _group(unfaulted, 3)
    _disable_faults(feature)
    np.testing.assert_allclose(values, _reference_value(feature, points), atol=1e-8, equal_nan=True)


# ---------------------------------------------------------------------------
# Original mask coordinates
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["onlap_domain", "stopped_fault"])
def test_region_masks_use_the_original_points(name, monkeypatch):
    """A region of a group gets the original points, never the restored points.

    A spy on every region of every group records the batch that the group mask
    gives to it. Calls inside ``apply_to_points`` (fault frames) are not recorded.
    The values must also equal the golden values. A wrong implementation that
    gives restored points to the region predicates fails both checks.
    """
    model, golden = build_model(name), _golden(name)
    points = golden["points"]
    batches: list[np.ndarray] = []
    depth = [0]
    original = FaultSegment.apply_to_points

    def tracked(self, pts, *args, **kwargs):
        depth[0] += 1
        try:
            return original(self, pts, *args, **kwargs)
        finally:
            depth[0] -= 1

    monkeypatch.setattr(FaultSegment, "apply_to_points", tracked)
    groups = [
        feature
        for feature in model.features
        if feature.name in model.feature_name_index and feature.name.startswith("group_")
    ]
    regions = {id(r): r for group in groups for r in group.regions if isinstance(r, BaseFeature)}
    assert regions
    for region in regions.values():
        region_value = region.evaluate_value

        def spy(pos, *args, _original=region_value, **kwargs):
            if depth[0] == 0:
                batches.append(np.array(pos))
            return _original(pos, *args, **kwargs)

        monkeypatch.setattr(region, "evaluate_value", spy)

    actual = {"points": points}
    for group in groups:
        actual[f"feature::{group.name}::value"] = group.evaluate_value(points)
    model_ids = model.evaluate_model(points)
    assert batches
    for batch in batches:
        np.testing.assert_array_equal(batch, points)
    _assert_golden(name, actual)
    _assert_golden(name, {"points": points, "model_ids": model_ids})
    # The restored points differ from the original points, so the check is real.
    fault_group = next(g for g in groups if g.faults)
    assert not np.allclose(_manual_restore(fault_group, points), points)


# ---------------------------------------------------------------------------
# Child-fault batches (recursive fault-frame batches, OP1)
# ---------------------------------------------------------------------------


def test_child_fault_frame_batches_are_counted_separately():
    """``cutting_faults``: group_0 uses [fault_1, fault_0]; the frame of fault_1 uses [fault_0].

    The frame of fault_1 has 3 coordinates and 10 steps. Distinct nested batches:
    the original points N (3 coordinates share it: 1 application), then the
    tetrahedron corners of the first step (1), then 9 further steps, each with
    one moved batch for the 3 values and one corner batch for the gradient (2 each).
    Nested: 1 + 1 + 9 * 2 = 20. Top level: fault_1 and fault_0 on the group
    batch = 2. Expected total 22. Legacy: 42 (40 nested + 2).
    """
    model, golden = build_model("cutting_faults"), _golden("cutting_faults")
    points = golden["points"]
    ids, counts = _counted(lambda: model.evaluate_model(points))
    _assert_golden("cutting_faults", {"points": points, "model_ids": ids})
    assert counts["fault_calls_total"] == 22
    assert counts["fault_calls_top_level"] == 2
    assert counts["fault_calls_total"] <= 42


# ---------------------------------------------------------------------------
# Opaque callables
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["coordinate", "feature"])
def test_opaque_region_callable_keeps_the_legacy_result(kind):
    """A plain function used as a region is opaque. Sharing is not assumed.

    Only the values are checked: they must equal the legacy recomputation, and no
    exception may happen.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 0)
    if kind == "coordinate":
        feature.add_region(lambda xyz: xyz[:, 0] < 600.0)
    else:
        # The region calls the evaluation of another faulted feature (group_3 is
        # eroded at some of the points). It is a callable, not a feature object.
        other = _group(model, 3)
        feature.add_region(lambda xyz: np.isnan(other.evaluate_value(xyz)))
    values = feature.evaluate_value(points)
    np.testing.assert_allclose(values, _reference_value(feature, points), atol=1e-8, equal_nan=True)
    assert np.isnan(values).any()
    assert np.isfinite(values).any()


# ---------------------------------------------------------------------------
# Exception cleanup
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("entry", ["value", "model"])
def test_exception_in_a_fault_application_resets_the_scope(entry, monkeypatch):
    """The second fault application raises. The next call must be complete and correct.

    A partial restoration of the first application must not stay in a scope. The
    next call therefore makes all 4 applications again and gives the golden values.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]

    def call():
        if entry == "value":
            return _group(model, 3).evaluate_value(points)
        return model.evaluate_model(points)

    original = FaultSegment.apply_to_points
    state = {"calls": 0}

    def failing(self, pts, *args, **kwargs):
        state["calls"] += 1
        if state["calls"] == 2:
            raise RuntimeError("injected failure")
        return original(self, pts, *args, **kwargs)

    with monkeypatch.context() as scoped:
        scoped.setattr(FaultSegment, "apply_to_points", failing)
        with pytest.raises(RuntimeError, match="injected failure"):
            call()
    result, counts = _counted(call)
    if entry == "value":
        _assert_golden("f4_u3", {"points": points, "feature::group_3::value": result})
    else:
        _assert_golden("f4_u3", {"points": points, "model_ids": result})
    assert counts["fault_calls_total"] == 4
    _assert_no_stale_scope(_group(model, 3), points, model)


def test_exception_after_a_complete_set_resets_the_scope(monkeypatch):
    """Set A is complete and the failure happens in set B (call 6 of 8).

    The model has 2 ordered fault sets, A and B, with 8 applications. A scope that is not
    reset in ``finally`` keeps set A, and a published partial set B is reused. Both give
    fewer than 8 applications in the next call. The values must equal an unscoped reference.
    """
    golden = _golden("f4_u3")
    points = golden["points"]
    reference = _reordered_model().evaluate_model(points)
    model = _reordered_model()
    original = FaultSegment.apply_to_points
    state = {"calls": 0}

    def failing(self, pts, *args, **kwargs):
        state["calls"] += 1
        if state["calls"] == 6:
            raise RuntimeError("injected failure")
        return original(self, pts, *args, **kwargs)

    with monkeypatch.context() as scoped:
        scoped.setattr(FaultSegment, "apply_to_points", failing)
        with pytest.raises(RuntimeError, match="injected failure"):
            model.evaluate_model(points)
    assert state["calls"] == 6
    ids, counts = _counted(lambda: model.evaluate_model(points))
    assert counts["fault_calls_total"] == 8
    np.testing.assert_array_equal(ids, reference)


# ---------------------------------------------------------------------------
# Standalone feature
# ---------------------------------------------------------------------------


def _detach_model(feature, seen: set[int] | None = None) -> None:
    """Remove the model from a feature and from all its region features."""
    seen = set() if seen is None else seen
    if id(feature) in seen or not isinstance(feature, BaseFeature):
        return
    seen.add(id(feature))
    feature._model = None
    for region in feature.regions:
        _detach_model(region, seen)


def test_standalone_feature_uses_a_temporary_scope():
    """A feature and its regions without a model share the restorations inside one call only.

    D3: a standalone feature uses a temporary scope. The regions have no model either, so
    their scope is the same temporary scope (a region with a model would be another scope).
    After the call no scope is active. Legacy: 16 applications.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points = golden["points"]
    feature = _group(model, 3)
    _detach_model(feature)
    assert feature.model is None
    assert all(region.model is None for region in feature.regions)
    values, counts = _counted(lambda: feature.evaluate_value(points))
    _assert_golden("f4_u3", {"points": points, "feature::group_3::value": values})
    assert counts["fault_calls_total"] == 4
    _assert_no_stale_scope(feature, points, None)


# ---------------------------------------------------------------------------
# Private API: evaluation_scope, EvaluationContext, restore
# ---------------------------------------------------------------------------


def test_scope_is_a_context_manager_that_gives_a_context():
    """``evaluation_scope(model=None)`` is a context manager that gives an EvaluationContext."""
    scope, context_type = _api("evaluation_scope"), _api("EvaluationContext")
    manager = scope()
    assert hasattr(manager, "__enter__") and hasattr(manager, "__exit__")
    with manager as context:
        assert isinstance(context, context_type)
    assert "model" in inspect.signature(scope).parameters
    assert inspect.signature(scope).parameters["model"].default is None


def test_nested_scope_reuses_the_context():
    """A nested scope for the same model gives the same context object."""
    scope = _api("evaluation_scope")
    model = build_model("f1_u0")
    with scope(model) as outer:
        with scope(model) as inner:
            assert inner is outer
        with scope(model) as again:
            assert again is outer
    with scope() as first, scope() as second:
        assert second is first


def test_scope_is_reset_after_exit_and_after_an_exception():
    """After the outer scope exits, also by an exception, a new scope is a new context."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, feature = golden["points"], _group(model, 3)
    with scope(model):
        restore(feature, points)
    # No stored entry stays: a restore after the exit applies the faults again.
    _assert_no_stale_scope(feature, points, model)
    with pytest.raises(RuntimeError, match="inside"), scope(model):
        restore(feature, points)
        raise RuntimeError("inside")
    _assert_no_stale_scope(feature, points, model)


def test_scope_of_another_model_is_separate():
    """A scope for model B inside a scope for model A is another context."""
    scope = _api("evaluation_scope")
    model_a, model_b = build_model("f1_u0"), build_model("f1_u0")
    with scope(model_a) as context_a:
        with scope(model_b) as context_b:
            assert context_b is not context_a
        with scope(model_a) as context_again:
            assert context_again is context_a


def test_restore_signature():
    """``restore(feature, points, reverse=False)`` is the cache entry of D3."""
    parameters = inspect.signature(_api("restore")).parameters
    assert list(parameters) == ["feature", "points", "reverse"]
    assert parameters["reverse"].default is False


def test_restore_twice_in_one_scope_applies_the_faults_once():
    """The second restore of equal points (also a copy) gives an equal array and no application."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, feature = golden["points"], _group(model, 3)
    expected = _manual_restore(feature, points)

    def run():
        with scope(model):
            first = restore(feature, points)
            second = restore(feature, points.copy())
        return first, second

    (first, second), counts = _counted(run)
    assert counts["fault_calls_total"] == 4
    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(first, expected)


def test_restored_array_is_read_only():
    """The shared array is read-only, so a consumer cannot change the shared value."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    with scope(model):
        restored = restore(_group(model, 3), golden["points"])
        assert restored.flags.writeable is False
        with pytest.raises(ValueError):
            restored[0, 0] = 0.0
    # The input array is not changed and stays writeable.
    assert golden["points"].flags.writeable is True


def test_restore_reverse_is_a_separate_batch_and_matches_the_legacy_result():
    """``reverse`` is part of the key: forward and reverse restorations are 2 x 4 applications."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, feature = golden["points"], _group(model, 3)

    def run():
        with scope(model):
            forward = restore(feature, points)
            backward = restore(feature, points, reverse=True)
            repeated = restore(feature, points, reverse=True)
        return forward, backward, repeated

    (forward, backward, repeated), counts = _counted(run)
    assert counts["fault_calls_total"] == 8
    np.testing.assert_array_equal(backward, repeated)
    np.testing.assert_array_equal(backward, _manual_restore(feature, points, reverse=True))
    assert not np.allclose(forward, backward)


def test_restore_with_disabled_faults_returns_the_points():
    """``faults_enabled = False`` gives equal points and no application."""
    scope, restore = _api("evaluation_scope"), _api("restore")
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, feature = golden["points"], _group(model, 3)
    feature.faults_enabled = False

    def run():
        with scope(model):
            return restore(feature, points)

    restored, counts = _counted(run)
    assert counts["fault_calls_total"] == 0
    np.testing.assert_array_equal(restored, points)


def test_public_calls_leave_no_scope_active():
    """After a top-level call, a new scope has no stored restoration (value and model value).

    The gradient call is only used to check that it leaves no scope active. No count is checked.
    """
    model, golden = build_model("f4_u3"), _golden("f4_u3")
    points, feature = golden["points"], _group(model, 3)
    feature.evaluate_value(points)
    _assert_no_stale_scope(feature, points, model)
    model.evaluate_model_gradient(points)
    _assert_no_stale_scope(feature, points, model)
