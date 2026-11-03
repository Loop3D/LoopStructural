"""Tests for the complete invalidation of the persistent restoration cache (D4, FR-009, FR-010).

A revision counter alone cannot find every change: public lists, dictionaries and arrays
are mutable. D4 therefore uses two guards. Both are private API, and this file fixes the
names for task T017.

``GeologicalModel.evaluation_revision``
    A read-only ``int`` property. It never decreases. A pure evaluation call does not change
    it. Every known mutation path of the objects of the model increases it: a setter or a
    method of a builder, an interpolator, a feature, a fault or the model.

``dependency_fingerprint(feature)`` in ``LoopStructural.modelling.features._evaluation_cache``
    A 32-byte SHA256 digest of all state that the restoration of the ordered faults of the
    feature depends on, or None when a reachable dependency is opaque. The function builds
    first (a dirty builder is built exactly once). It never changes the model. It ends on a
    graph with cycles. It hashes arrays in chunks and does not copy a whole array.

Test groups
-----------
* Mutation tests use the golden models. These models use stored fields, so no test of this
  group starts a solve. Each mutation must change the fingerprint, must change the model
  result (the guard), and must give a fresh result when the cache is on.
* Builder tests use a TINY real model. It has one fault, one group and very few data points,
  so a real build takes a few hundred milliseconds.
* Opaque, cycle and memory tests check the rules for the cases where no digest is possible.

The tests in the group "mutation changes the model result" use only the present API. They
pass before T017. All other tests fail while the new API is absent.
"""

from __future__ import annotations

import importlib
import threading
import tracemalloc
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from loop_interpolation._discrete_interpolator import DiscreteInterpolator

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation
from LoopStructural import GeologicalModel
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn
from LoopStructural.modelling.features._base_geological_feature import BaseFeature
from LoopStructural.modelling.features.builders._geological_feature_builder import (
    GeologicalFeatureBuilder,
)
from LoopStructural.utils import PositiveRegion
from tests.fault_evaluation.models import build_model

MODULE_NAME = "LoopStructural.modelling.features._evaluation_cache"
DATA = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _api(name: str):
    """Return an object of the cache module, or fail with a clear reason when it is absent."""
    try:
        module = importlib.import_module(MODULE_NAME)
    except ImportError as error:
        pytest.fail(f"{MODULE_NAME} is required by D4 (import failed: {error})")
    value = getattr(module, name, None)
    assert value is not None, f"{MODULE_NAME}.{name} is required by D4"
    return value


def _revision(model: GeologicalModel) -> int:
    """Return ``model.evaluation_revision``, or fail with a clear reason when it is absent."""
    value = getattr(model, "evaluation_revision", None)
    assert value is not None, "GeologicalModel.evaluation_revision is required by D4"
    assert type(value) is int and value >= 0, f"revision must be a non-negative int, got {value!r}"
    return value


def _raw_fingerprint(feature) -> bytes | None:
    """Call ``dependency_fingerprint`` and return its result (bytes or None)."""
    result = _api("dependency_fingerprint")(feature)
    assert result is None or (
        isinstance(result, bytes) and len(result) == 32
    ), f"dependency_fingerprint must return None or a 32-byte digest, got {result!r}"
    return result


def _fingerprint(feature) -> bytes:
    """Return the digest of a feature that must be supported (not opaque)."""
    result = _raw_fingerprint(feature)
    assert result is not None, f"feature {feature.name!r} must not be opaque"
    return result


def _points(name: str) -> np.ndarray:
    """Return the stored golden points of a model."""
    with np.load(DATA / f"{name}.npz", allow_pickle=False) as archive:
        return archive["points"].copy()


def _evaluate_disabled(model: GeologicalModel, call: Callable[[], np.ndarray]):
    """Run ``call`` with the cache off. Return the result and the fault call counters."""
    cache = model.evaluation_cache
    previous = cache.enabled
    cache.enabled = False
    try:
        with EvaluationInstrumentation() as instrumentation:
            result = call()
        return result, instrumentation.snapshot()
    finally:
        cache.enabled = previous


def _fault(model: GeologicalModel, index: int):
    return model.get_feature_by_name(f"fault_{index}")


def _group(model: GeologicalModel, index: int):
    return model.get_feature_by_name(f"group_{index}")


def _assign(obj, name: str, value) -> Callable[[], None]:
    """Set an attribute. Return the function that sets the old value back."""
    old = getattr(obj, name)
    setattr(obj, name, value)
    return lambda: setattr(obj, name, old)


# ---------------------------------------------------------------------------
# Mutations of the golden models. Each function returns its undo function.
# ---------------------------------------------------------------------------


def _m_faults_reorder(model):
    group = _group(model, 3)
    return _assign(group, "faults", list(reversed(group.faults)))


def _m_faults_subset(model):
    group = _group(model, 3)
    return _assign(group, "faults", list(group.faults[:2]))


def _m_faults_inplace_reverse(model):
    faults = _group(model, 3).faults
    faults.reverse()
    return faults.reverse


def _m_faults_inplace_pop(model):
    faults = _group(model, 3).faults
    removed = faults.pop()
    return lambda: faults.append(removed)


def _m_faults_enabled(model):
    return _assign(_group(model, 3), "faults_enabled", False)


def _m_displacement_attribute(model):
    return _assign(_fault(model, 0), "displacement", 20.0)


def _m_displacement_method(model):
    fault = _fault(model, 1)
    old = fault.displacement
    fault.set_displacement(20.0, scale=False)
    return lambda: fault.set_displacement(old, scale=False)


def _m_steps(model):
    return _assign(_fault(model, 2), "steps", 2)


def _m_faultfunction_setter(model):
    return _assign(_fault(model, 3), "faultfunction", "BaseFault3D")


def _m_faultfunction_inplace(model):
    # The default function is one object that all faults share. The undo is mandatory.
    function = _fault(model, 0).faultfunction
    old = function.scale
    function.scale = 0.25

    def undo():
        function.scale = old

    return undo


def _m_region_add(model):
    group = _group(model, 3)
    region = PositiveRegion(
        _fault(model, 0),
        vector=np.array([1.0, 0.0, 0.0]),
        point=np.array([500.0, 500.0, -500.0]),
    )
    group.add_region(region)
    return lambda: group.regions.remove(region)


def _m_region_remove_inplace(model):
    regions = _group(model, 3).regions
    removed = regions.pop(2)
    return lambda: regions.insert(2, removed)


def _m_region_threshold(model):
    return _assign(_group(model, 3).regions[2], "value", 50.0)


def _m_region_sign(model):
    return _assign(_group(model, 3).regions[2], "sign", False)


def _c_inplace(fault_index: int) -> Callable:
    def mutate(model):
        coefficients = _fault(model, fault_index).features[0].interpolator.c
        old = coefficients.copy()
        coefficients += 0.3 * np.nanstd(coefficients)

        def undo():
            coefficients[:] = old

        return undo

    return mutate


def _m_c_rebind(model):
    interpolator = _fault(model, 0).features[0].interpolator
    return _assign(interpolator, "c", interpolator.c * 1.2)


def _support(model):
    return _fault(model, 0).features[0].interpolator.support


def _m_support_origin(model):
    origin = _support(model).origin
    old = origin.copy()
    origin += 30.0

    def undo():
        origin[:] = old

    return undo


def _m_support_step(model):
    step = _support(model).step_vector
    old = step.copy()
    step *= 1.05

    def undo():
        step[:] = old

    return undo


def _m_support_rotation(model):
    support = _support(model)
    old = support.rotation_xy.copy()
    support.rotation_xy = 15.0
    return lambda: setattr(support, "rotation_xy", old)


@dataclass(frozen=True)
class Mutation:
    """One change of the golden model.

    ``revision`` is True for a setter or a method path: the revision must increase. It is
    False for an in-place edit: only the fingerprint must change. ``effect`` is False for
    a change that does not change the model result today (bookkeeping): the guard is then
    skipped.
    """

    key: str
    model: str
    target: str
    apply: Callable
    revision: bool
    effect: bool = True


MUTATIONS = [
    Mutation("faults_reorder", "f4_u3", "group_3", _m_faults_reorder, True),
    Mutation("faults_subset", "f4_u3", "group_3", _m_faults_subset, True),
    Mutation("faults_inplace_reverse", "f4_u3", "group_3", _m_faults_inplace_reverse, False),
    Mutation("faults_inplace_pop", "f4_u3", "group_3", _m_faults_inplace_pop, False),
    Mutation("faults_enabled", "f4_u3", "group_3", _m_faults_enabled, True),
    Mutation("displacement_attribute", "f4_u3", "group_3", _m_displacement_attribute, True),
    Mutation("displacement_method", "f4_u3", "group_3", _m_displacement_method, True),
    Mutation("steps", "f4_u3", "group_3", _m_steps, True),
    Mutation("faultfunction_setter", "f4_u3", "group_3", _m_faultfunction_setter, True),
    Mutation("faultfunction_inplace", "f4_u3", "group_3", _m_faultfunction_inplace, False),
    Mutation("region_add", "f4_u3", "group_3", _m_region_add, True),
    Mutation("region_remove_inplace", "f4_u3", "group_3", _m_region_remove_inplace, False),
    Mutation("region_threshold", "f4_u3", "group_3", _m_region_threshold, False),
    Mutation("region_sign", "f4_u3", "group_3", _m_region_sign, False),
    Mutation("frame_c_inplace", "f4_u3", "group_3", _c_inplace(0), False),
    Mutation("frame_c_rebind", "f4_u3", "group_3", _m_c_rebind, False),
    Mutation("frame_c_inplace_child", "cutting_faults", "group_0", _c_inplace(1), False),
    Mutation("frame_c_inplace_parent", "cutting_faults", "group_0", _c_inplace(0), False),
    Mutation("support_origin_inplace", "f4_u3", "group_3", _m_support_origin, False),
    Mutation("support_step_inplace", "f4_u3", "group_3", _m_support_step, False),
    # StructuredGrid evaluation does not use rotation_xy today: bookkeeping for the guard.
    Mutation("support_rotation", "f4_u3", "group_3", _m_support_rotation, False, effect=False),
]
MUTATION_IDS = [mutation.key for mutation in MUTATIONS]


# ---------------------------------------------------------------------------
# Group 1: a mutation changes the model result (present API only, passes before T017)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mutation", MUTATIONS, ids=MUTATION_IDS)
def test_mutated_model_gives_a_fresh_changed_result(mutation):
    """After a mutation the cached evaluation equals the cache-disabled evaluation.

    The first evaluation warms the cache. After the mutation the cached evaluation must
    apply as many faults as the cache-disabled one, so no stale lookup is used. The guard
    shows that the mutation really changes the output: a mutation without effect would
    let every other test pass for a wrong reason.
    """
    model = build_model(mutation.model)
    feature = model.get_feature_by_name(mutation.target)
    points = _points(mutation.model)
    before = feature.evaluate_value(points)
    model.evaluate_model(points)
    undo = mutation.apply(model)
    try:
        with EvaluationInstrumentation() as instrumentation:
            cached = feature.evaluate_value(points)
        cached_counts = instrumentation.snapshot()
        cached_ids = model.evaluate_model(points)
        disabled, disabled_counts = _evaluate_disabled(
            model, lambda: feature.evaluate_value(points)
        )
        disabled_ids, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
        np.testing.assert_array_equal(cached, disabled)
        np.testing.assert_array_equal(cached_ids, disabled_ids)
        assert cached_counts["fault_calls_total"] == disabled_counts["fault_calls_total"]
        if mutation.effect:
            assert not np.array_equal(disabled, before, equal_nan=True), "no effect"
    finally:
        undo()


# ---------------------------------------------------------------------------
# Group 2: fingerprint and revision for the golden models
# ---------------------------------------------------------------------------


def test_api_has_the_names_fixed_for_t017():
    """``dependency_fingerprint`` takes one feature. The revision is a read-only property."""
    function = _api("dependency_fingerprint")
    assert callable(function)
    model = build_model("f1_u0")
    revision = _revision(model)
    prop = GeologicalModel.__dict__.get("evaluation_revision")
    assert isinstance(prop, property) and prop.fset is None
    with pytest.raises(AttributeError):
        model.evaluation_revision = revision + 1
    assert _revision(model) == revision
    digest = _fingerprint(_group(model, 0))
    assert isinstance(digest, bytes) and len(digest) == 32


@pytest.mark.parametrize("mutation", MUTATIONS, ids=MUTATION_IDS)
def test_mutation_changes_the_fingerprint(mutation):
    """Each mutation changes the digest. A setter or method path also increases the revision.

    A second call on the changed model gives the same digest, and the call does not change
    the revision. When the change is undone, the original digest returns, because the
    output depends on the full current state. The revision never goes back.
    """
    model = build_model(mutation.model)
    feature = model.get_feature_by_name(mutation.target)
    before, revision = _fingerprint(feature), _revision(model)
    undo = mutation.apply(model)
    try:
        after, revised = _fingerprint(feature), _revision(model)
        assert after != before, "the fingerprint did not change"
        if mutation.revision:
            assert revised > revision, "the revision did not increase"
        else:
            assert revised >= revision
        assert _fingerprint(feature) == after
        assert _revision(model) == revised
    finally:
        undo()
    assert _fingerprint(feature) == before, "the restored state must give the first digest"
    assert _revision(model) >= revised


def _noop_evaluate_value(model):
    _group(model, 3).evaluate_value(_points("f4_u3"))


def _noop_evaluate_model(model):
    model.evaluate_model(_points("f4_u3"))


def _noop_evaluate_model_gradient(model):
    model.evaluate_model_gradient(_points("f4_u3"))


def _noop_evaluate_feature_value(model):
    model.evaluate_feature_value("group_3", _points("f4_u3"))


def _noop_fingerprint_call(model):
    _fingerprint(_group(model, 3))
    _fingerprint(_fault(model, 0))


def _noop_same_displacement(model):
    fault = _fault(model, 0)
    fault.displacement = fault.displacement


def _noop_same_steps(model):
    fault = _fault(model, 0)
    fault.steps = fault.steps


def _noop_same_faults_enabled(model):
    group = _group(model, 3)
    group.faults_enabled = group.faults_enabled


def _noop_same_faults(model):
    group = _group(model, 3)
    group.faults = list(group.faults)


def _noop_same_faultfunction(model):
    fault = _fault(model, 0)
    fault.faultfunction = fault.faultfunction


def _noop_same_build_arguments(model):
    builder = _group(model, 3).builder
    builder.update_build_arguments(dict(builder.build_arguments))


NOOPS = [
    _noop_evaluate_value,
    _noop_evaluate_model,
    _noop_evaluate_model_gradient,
    _noop_evaluate_feature_value,
    _noop_fingerprint_call,
    _noop_same_build_arguments,
]

# A setter may always increase the revision, also for an equal value. The plan (D4) does not
# forbid this. For these cases only the fingerprint is checked.
SAME_VALUE_ASSIGNMENTS = [
    _noop_same_displacement,
    _noop_same_steps,
    _noop_same_faults_enabled,
    _noop_same_faults,
    _noop_same_faultfunction,
]


@pytest.mark.parametrize("noop", NOOPS, ids=[noop.__name__[6:] for noop in NOOPS])
def test_noop_keeps_the_fingerprint_and_the_revision(noop):
    """A pure evaluation and a digest call change neither the fingerprint nor the revision."""
    model = build_model("f4_u3")
    features = [_group(model, 3), _group(model, 0), _fault(model, 0)]
    before = [_fingerprint(feature) for feature in features]
    revision = _revision(model)
    noop(model)
    assert [_fingerprint(feature) for feature in features] == before
    assert _revision(model) == revision


@pytest.mark.parametrize(
    "assign", SAME_VALUE_ASSIGNMENTS, ids=[assign.__name__[6:] for assign in SAME_VALUE_ASSIGNMENTS]
)
def test_same_value_assignment_keeps_the_fingerprint(assign):
    """An assignment of an equal value does not change the fingerprint.

    The revision is not checked. A setter may increase it, also for an equal value.
    """
    model = build_model("f4_u3")
    features = [_group(model, 3), _group(model, 0), _fault(model, 0)]
    before = [_fingerprint(feature) for feature in features]
    assign(model)
    assert [_fingerprint(feature) for feature in features] == before


def test_revision_is_a_read_only_int_that_never_decreases():
    """The property gives an int, rejects assignment, and the value only grows."""
    model = build_model("f4_u3")
    seen = [_revision(model)]
    with pytest.raises(AttributeError):
        model.evaluation_revision = 0
    group, fault = _group(model, 3), _fault(model, 0)
    for step in (
        lambda: setattr(fault, "displacement", 10.0),
        lambda: group.evaluate_value(_points("f4_u3")),
        lambda: setattr(fault, "displacement", 40.0),
        lambda: setattr(group, "faults_enabled", False),
        lambda: setattr(group, "faults_enabled", True),
        lambda: setattr(fault, "steps", 3),
        lambda: setattr(fault, "steps", 10),
    ):
        step()
        seen.append(_revision(model))
    assert seen == sorted(seen), f"the revision decreased: {seen}"
    assert seen[-1] > seen[0]


# ---------------------------------------------------------------------------
# Unrelated features and dependency frames
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("change", ["displacement", "steps", "frame_c"])
def test_unrelated_feature_fingerprint_stays_equal(change):
    """``stopped_fault``: group_0 has no fault. A change of fault_0 changes group_1 only."""
    model = build_model("stopped_fault")
    group_0, group_1, fault = _group(model, 0), _group(model, 1), _fault(model, 0)
    assert group_0.faults == [] and group_1.faults == [fault]
    before_0, before_1 = _fingerprint(group_0), _fingerprint(group_1)
    if change == "displacement":
        fault.displacement = 20.0
    elif change == "steps":
        fault.steps = 3
    else:
        fault.features[0].interpolator.c[:] += 1.0
    assert _fingerprint(group_0) == before_0
    assert _fingerprint(group_1) != before_1


def test_changing_one_fault_leaves_the_other_independent_faults_equal():
    """``f4_u3``: the 4 faults do not reach each other. A change of fault_3 keeps the others."""
    model = build_model("f4_u3")
    faults = [_fault(model, index) for index in range(4)]
    before = [_fingerprint(fault) for fault in faults]
    faults[3].displacement = 15.0
    after = [_fingerprint(fault) for fault in faults]
    assert after[:3] == before[:3]


FRAME_CASES = [
    # model, mutated fault, features that must change, features that must stay equal
    ("cutting_faults", 0, ["fault_1", "group_0"], ["fault_0"]),
    ("cutting_faults", 1, ["group_0"], ["fault_0"]),
    ("f4_u3", 2, ["group_0", "group_1", "group_2", "group_3"], ["fault_0", "fault_1", "fault_3"]),
]


@pytest.mark.parametrize("name,index,changed,unchanged", FRAME_CASES)
def test_frame_coefficient_changes_every_feature_that_reaches_the_fault(
    name, index, changed, unchanged
):
    """A change of a coefficient of a fault frame changes all features that reach the fault.

    ``cutting_faults``: group_0 uses [fault_1, fault_0] and the frame of fault_1 is faulted
    by fault_0. The guard shows that the output of group_0 really changes.
    """
    model = build_model(name)
    points = _points(name)
    before = {key: _fingerprint(model.get_feature_by_name(key)) for key in changed + unchanged}
    values, _ = _evaluate_disabled(
        model, lambda: model.get_feature_by_name("group_0").evaluate_value(points)
    )
    revision = _revision(model)
    coefficients = _fault(model, index).features[0].interpolator.c
    coefficients += 0.3 * np.nanstd(coefficients)
    for key in changed:
        assert _fingerprint(model.get_feature_by_name(key)) != before[key], key
    for key in unchanged:
        assert _fingerprint(model.get_feature_by_name(key)) == before[key], key
    assert _revision(model) >= revision
    changed_values, _ = _evaluate_disabled(
        model, lambda: model.get_feature_by_name("group_0").evaluate_value(points)
    )
    assert not np.array_equal(values, changed_values, equal_nan=True), "no effect"


@pytest.mark.parametrize("name", ["f4_u3", "cutting_faults", "stopped_fault", "f1_u0"])
def test_golden_features_are_not_opaque(name):
    """All objects of these models are known types. No feature gives None."""
    model = build_model(name)
    for feature in model.features:
        assert _raw_fingerprint(feature) is not None, f"{name}:{feature.name} is opaque"


# ---------------------------------------------------------------------------
# Opaque dependencies: no digest, so no persistent reuse
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("where", ["feature", "reachable_fault"])
@pytest.mark.parametrize("kind", ["lambda", "closure"])
def test_opaque_region_predicate_gives_none(where, kind):
    """A plain function used as a region has unknown state. The result is None.

    The function may read a closure or a global that changes without notice, so an object
    id never proves that it is unchanged. Without the function, the digest returns.
    """
    model = build_model("cutting_faults")
    group, fault = _group(model, 0), _fault(model, 0)
    before, revision = _fingerprint(group), _revision(model)
    state = {"limit": 600.0}
    if kind == "lambda":

        def predicate(xyz):
            return xyz[:, 0] < 600.0

    else:

        def predicate(xyz):
            return xyz[:, 0] < state["limit"]

    owner = group if where == "feature" else fault
    owner.regions.append(predicate)
    revision = _revision(model)
    assert _raw_fingerprint(group) is None
    state["limit"] = 100.0
    assert _raw_fingerprint(group) is None
    assert _revision(model) == revision, "dependency_fingerprint must not change the model"
    owner.regions.remove(predicate)
    assert _fingerprint(group) == before


class _CustomInterpolator:
    """An interpolator that is not a DiscreteInterpolator. Its state is unknown."""

    up_to_date = True

    def __init__(self, real):
        self.c = real.c
        self.support = real.support
        self.data = {}
        self._real = real

    def evaluate_value(self, xyz):
        return self._real.evaluate_value(xyz)

    def evaluate_gradient(self, xyz):
        return self._real.evaluate_gradient(xyz)


class _UnknownSupport:
    """A support of a type that the fingerprint does not know."""


@pytest.mark.parametrize("kind", ["custom_interpolator", "unsupported_support"])
def test_unsupported_interpolator_or_support_gives_none(kind):
    """A custom interpolator and an unknown support are opaque. The digest comes back later."""
    model = build_model("cutting_faults")
    group = _group(model, 0)
    coordinate = _fault(model, 0).features[0]
    before = _fingerprint(group)
    old_feature_interpolator = coordinate.interpolator
    old_builder_interpolator = coordinate.builder._interpolator
    old_support = old_feature_interpolator.support
    try:
        if kind == "custom_interpolator":
            fake = _CustomInterpolator(old_feature_interpolator)
            coordinate.interpolator = fake
            coordinate.builder._interpolator = fake
        else:
            old_feature_interpolator.support = _UnknownSupport()
        assert _raw_fingerprint(group) is None
    finally:
        coordinate.interpolator = old_feature_interpolator
        coordinate.builder._interpolator = old_builder_interpolator
        old_feature_interpolator.support = old_support
    assert _fingerprint(group) == before


# ---------------------------------------------------------------------------
# Cycles
# ---------------------------------------------------------------------------


class _Node(BaseFeature):
    """A small feature without builder. Tests link these nodes to make cycles."""

    def evaluate_value(self, pos, ignore_regions=False):
        return np.zeros(len(pos))

    def evaluate_gradient(self, pos, ignore_regions=False):
        return np.zeros((len(pos), 3))

    def get_data(self, value_map=None):
        return []

    def copy(self, name=None):
        return _Node(name or self.name)


def _cycle(kind: str) -> _Node:
    """Return a node that is part of a cycle of the chosen kind."""
    a, b, c = _Node("a"), _Node("b"), _Node("c")
    if kind == "faults":
        a._faults, b._faults, c._faults = [b], [c], [a]
    elif kind == "self_fault":
        a._faults = [a]
    elif kind == "frame":
        a.features, b.features, c.features = [b], [c], [a]
    elif kind == "parent":
        a.parent, b.parent, c.parent = b, c, a
    elif kind == "self_parent":
        a.parent = a
    elif kind == "regions":
        a.regions, b.regions = [b], [a]
    else:
        raise ValueError(kind)
    return a


def _run_with_timeout(call: Callable, seconds: float = 30.0):
    """Run ``call`` in a thread. Fail when it does not end, re-raise its error."""
    outcome: dict = {}

    def target():
        try:
            outcome["value"] = call()
        except BaseException as error:  # noqa: BLE001 - the error is re-raised below
            outcome["error"] = error

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(seconds)
    assert not thread.is_alive(), "dependency_fingerprint did not end (cycle?)"
    if "error" in outcome:
        raise outcome["error"]
    return outcome["value"]


@pytest.mark.parametrize(
    "kind", ["faults", "self_fault", "frame", "parent", "self_parent", "regions"]
)
def test_cycle_in_the_feature_graph_ends(kind):
    """A graph with a cycle ends. The result is a digest or None, never an endless loop."""
    function = _api("dependency_fingerprint")
    node = _cycle(kind)
    result = _run_with_timeout(lambda: function(node))
    assert result is None or (isinstance(result, bytes) and len(result) == 32)


def test_real_self_referential_region_ends_with_a_digest():
    """In ``stopped_fault`` the inverse unconformity of group_0 has group_0 as its parent."""
    model = build_model("stopped_fault")
    group = _group(model, 0)
    assert any(getattr(region, "parent", None) is group for region in group.regions)
    digest = _run_with_timeout(lambda: _raw_fingerprint(group))
    assert digest is not None


# ---------------------------------------------------------------------------
# Memory
# ---------------------------------------------------------------------------


def test_fingerprint_does_not_copy_a_large_coefficient_array():
    """The hash reads the array in chunks. The peak of new memory stays far below its size.

    The array has 24 MB. A copy (for example ``tobytes``) would need 24 MB. The limit of
    6 MB leaves room for the chunk buffer, the digest state and Python objects.
    """
    model = build_model("f1_u0")
    group = _group(model, 0)
    interpolator = _fault(model, 0).features[0].interpolator
    interpolator.c = np.linspace(0.0, 1.0, 3_000_000)
    size = interpolator.c.nbytes
    _fingerprint(group)  # warm call: imports and lazy tables are not counted
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        _fingerprint(group)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 6 * 1024 * 1024 < size, f"peak {peak} bytes for an array of {size} bytes"


# ---------------------------------------------------------------------------
# Tiny real model: builders, data, build arguments, solves
# ---------------------------------------------------------------------------

TINY_NELEMENTS = 300


def _tiny_data() -> pd.DataFrame:
    """Return few data points: 1 fault (3 values, 1 normal) and 2 groups (value levels)."""
    rows = []
    for y in (100.0, 500.0, 900.0):
        rows.append({"X": 500.0, "Y": y, "Z": -500.0, "val": 0.0, "feature_name": "fault_0"})
    rows.append(
        {
            "X": 500.0,
            "Y": 500.0,
            "Z": -500.0,
            "gx": 1.0,
            "gy": 0.0,
            "gz": 0.0,
            "feature_name": "fault_0",
        }
    )
    rng = np.random.default_rng(3)
    for name, shift in (("group_0", -900.0), ("group_1", -300.0)):
        for value in (0.0, 50.0, 100.0):
            for x, y in rng.random((4, 2)) * 1000.0:
                rows.append(
                    {
                        "X": x,
                        "Y": y,
                        "Z": shift + value + 0.05 * x,
                        "val": value,
                        "feature_name": name,
                    }
                )
        rows.append(
            {
                "X": 500.0,
                "Y": 500.0,
                "Z": shift,
                "gx": 0.0,
                "gy": 0.0,
                "gz": 1.0,
                "feature_name": name,
            }
        )
    return pd.DataFrame(rows)


def _tiny_model() -> GeologicalModel:
    """Return a model with one fault and one group. Nothing is built yet (all is dirty)."""
    model = GeologicalModel(np.array([0.0, 0.0, -1000.0]), np.array([1000.0, 1000.0, 0.0]))
    model.data = _tiny_data()
    model.create_and_add_fault("fault_0", 40.0, nelements=TINY_NELEMENTS)
    model.create_and_add_foliation("group_0", nelements=TINY_NELEMENTS)
    return model


def _tiny_points() -> np.ndarray:
    return np.random.default_rng(1).random((60, 3)) * np.array([1000.0, 1000.0, -1000.0])


def _built_tiny_model():
    """Return a built tiny model, its group, its points and the output before any change."""
    model = _tiny_model()
    group, points = _group(model, 0), _tiny_points()
    values, _ = _evaluate_disabled(model, lambda: group.evaluate_value(points))
    return model, group, points, values


def _install_counters(monkeypatch):
    """Count the builds per builder and the solves, and call the real functions."""
    builds: dict[int, int] = {}
    solves = [0]
    original_build = GeologicalFeatureBuilder.build
    original_solve = DiscreteInterpolator.solve_system

    def counting_build(self, *args, **kwargs):
        builds[id(self)] = builds.get(id(self), 0) + 1
        return original_build(self, *args, **kwargs)

    def counting_solve(self, *args, **kwargs):
        solves[0] += 1
        return original_solve(self, *args, **kwargs)

    monkeypatch.setattr(GeologicalFeatureBuilder, "build", counting_build)
    monkeypatch.setattr(DiscreteInterpolator, "solve_system", counting_solve)
    return builds, solves


@pytest.mark.parametrize("case", ["all_dirty", "feature_builder_only", "fault_frame_builder_only"])
def test_dirty_lazy_build_happens_exactly_once(case, monkeypatch):
    """The first digest builds each dirty builder once. A second digest builds nothing.

    The build runs before the faults are read: a fresh feature has no fault list until its
    builder has run. The revision grows with the build and stays equal afterwards.
    """
    builds, solves = _install_counters(monkeypatch)
    model = _tiny_model()
    group, fault, points = _group(model, 0), _fault(model, 0), _tiny_points()
    frame_builders = list(fault.builder.builders)
    if case == "all_dirty":
        assert group.faults == []
        expected = {id(group.builder)} | {id(builder) for builder in frame_builders}
    else:
        group.evaluate_value(points)
        if case == "feature_builder_only":
            group.builder.set_not_up_to_date(None)
            expected = {id(group.builder)}
        else:
            frame_builders[1].set_not_up_to_date(None)
            expected = {id(frame_builders[1])}
    builds.clear()
    solves[0] = 0
    revision = _revision(model)
    first = _fingerprint(group)
    assert set(builds) == expected and set(builds.values()) == {1}, builds
    assert solves[0] == len(expected)
    assert group.faults == [fault], "the digest must be made after the build"
    assert _revision(model) > revision
    revision = _revision(model)
    second = _fingerprint(group)
    assert second == first
    assert sum(builds.values()) == len(expected) and solves[0] == len(expected)
    assert _revision(model) == revision
    _evaluate_disabled(model, lambda: group.evaluate_value(points))
    assert sum(builds.values()) == len(expected), "an evaluation after the digest must not build"


def _data_value_in_place(builder):
    builder.data.loc[builder.data["val"] == 50.0, "val"] = 60.0


def _data_x_in_place(builder):
    builder.data.loc[builder.data.index[:3], "X"] += 25.0


def _data_frame_replaced(builder):
    builder.data = builder.data.assign(Z=builder.data["Z"] + 30.0)


def _data_add_from_frame(builder):
    changed = builder.data.copy()
    changed.loc[changed["val"] == 50.0, "val"] = 60.0
    builder.add_data_from_data_frame(changed)


DATA_CASES = [
    ("value_in_place", _data_value_in_place, False),
    ("x_in_place", _data_x_in_place, False),
    ("frame_replaced", _data_frame_replaced, False),
    ("add_data_from_data_frame", _data_add_from_frame, True),
]


@pytest.mark.parametrize(
    "apply,revision",
    [(apply, revision) for _, apply, revision in DATA_CASES],
    ids=[key for key, _, _ in DATA_CASES],
)
def test_builder_data_change_changes_the_fingerprint(apply, revision):
    """A change of the data of a builder changes the digest, before any rebuild.

    An in-place edit and a direct assignment do not mark the builder dirty. The digest
    must still change. The method ``add_data_from_data_frame`` also increases the
    revision. The guard marks the builder dirty and shows that the new data change the
    output.
    """
    model, group, points, before_values = _built_tiny_model()
    before, rev = _fingerprint(group), _revision(model)
    apply(group.builder)
    assert _fingerprint(group) != before, "the fingerprint did not change"
    if revision:
        assert _revision(model) > rev
    group.builder.set_not_up_to_date(None)
    after_values, _ = _evaluate_disabled(model, lambda: group.evaluate_value(points))
    assert not np.array_equal(after_values, before_values, equal_nan=True), "no effect"


@pytest.mark.parametrize("key,value", [("regularisation", 0.05), ("cpw", 5.0), ("nelements", 150)])
def test_build_argument_change_changes_the_fingerprint(key, value):
    """``update_build_arguments`` increases the revision, marks dirty and changes the digest.

    The digest call builds once with the new argument. The result differs from the old one.
    """
    model, group, points, before_values = _built_tiny_model()
    before, rev = _fingerprint(group), _revision(model)
    group.builder.update_build_arguments({key: value})
    assert _revision(model) > rev
    assert _fingerprint(group) != before
    if key != "nelements":  # the builder stores the real element count for nelements
        assert group.builder.build_arguments[key] == value
    after_values, _ = _evaluate_disabled(model, lambda: group.evaluate_value(points))
    assert not np.array_equal(after_values, before_values, equal_nan=True), "no effect"


def _dirty(model, group):
    group.builder.set_not_up_to_date(None)


def _make_interpolator_dirty(model, group):
    group.interpolator.up_to_date = False


def _other_data(model) -> pd.DataFrame:
    changed = model.data.copy()
    changed.loc[changed["val"] == 50.0, "val"] = 55.0
    return changed


# (key, prepare, action). The revision is read after prepare, before action.
REVISION_PATHS = [
    ("set_not_up_to_date", None, lambda m, g: g.builder.set_not_up_to_date(None)),
    (
        "update_build_arguments",
        None,
        lambda m, g: g.builder.update_build_arguments({"regularisation": 0.5}),
    ),
    ("builder_build", None, lambda m, g: g.builder.build(**g.builder.build_arguments)),
    ("builder_update", _dirty, lambda m, g: g.builder.update()),
    ("builder_up_to_date", _dirty, lambda m, g: g.builder.up_to_date()),
    ("interpolator_solve_system", None, lambda m, g: g.interpolator.solve_system(solver="cg")),
    ("interpolator_update", _make_interpolator_dirty, lambda m, g: g.interpolator.update()),
    ("model_update", _dirty, lambda m, g: m.update(progressbar=False)),
    ("model_data", None, lambda m, g: setattr(m, "data", _other_data(m))),
    (
        "model_stratigraphic_column",
        None,
        lambda m, g: setattr(m, "stratigraphic_column", StratigraphicColumn()),
    ),
    (
        "model_add_feature",
        None,
        lambda m, g: m.create_and_add_foliation("group_1", nelements=TINY_NELEMENTS),
    ),
]


@pytest.mark.parametrize(
    "prepare,action",
    [(prepare, action) for _, prepare, action in REVISION_PATHS],
    ids=[key for key, _, _ in REVISION_PATHS],
)
def test_build_solve_and_model_paths_increase_the_revision(prepare, action):
    """Each known build, solve or model-level path increases the revision.

    These paths do not need an output guard: the change of data of the model, of the
    column and of the relations is bookkeeping for the stored restorations. A rebuild
    with equal data gives the same field up to the solver tolerance.
    """
    model, group, _points_unused, _values = _built_tiny_model()
    if prepare is not None:
        prepare(model, group)
    revision = _revision(model)
    action(model, group)
    assert _revision(model) > revision


def test_set_not_up_to_date_and_rebuild_with_equal_data_keep_the_result():
    """Bookkeeping case: dirty plus rebuild with equal data changes the result only by the solver tolerance."""
    model, group, points, before = _built_tiny_model()
    revision = _revision(model)
    group.builder.set_not_up_to_date(None)
    assert _revision(model) > revision
    after, _ = _evaluate_disabled(model, lambda: group.evaluate_value(points))
    np.testing.assert_allclose(after, before, rtol=0.0, atol=1e-2)


def test_revision_never_decreases_over_a_build_sequence():
    """A sequence of builds, evaluations and changes gives a non-decreasing revision."""
    model = _tiny_model()
    group, points = _group(model, 0), _tiny_points()
    seen = [_revision(model)]
    steps = [
        lambda: group.evaluate_value(points),
        lambda: group.evaluate_value(points),
        lambda: group.builder.update_build_arguments({"cpw": 2.0}),
        lambda: group.evaluate_value(points),
        lambda: _fingerprint(group),
        lambda: setattr(_fault(model, 0), "displacement", 10.0),
        lambda: group.evaluate_value(points),
    ]
    for step in steps:
        step()
        seen.append(_revision(model))
    assert seen == sorted(seen), f"the revision decreased: {seen}"
    # Two evaluations of an up-to-date model do not change it.
    assert seen[1] == seen[2]
