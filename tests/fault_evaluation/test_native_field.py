"""Native values, gradients and model ids agree with the Python path (spec 007, story 3).

Covers FR-004, FR-005, FR-005a and FR-005b. The native result is ``capture`` with the
defaults. The reference is ``capture`` with ``LoopStructuralConfig.native_faults`` off.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

pytest.importorskip("loop_model_core")

from LoopStructural import LoopStructuralConfig
from LoopStructural.modelling.features import _native_faults
from tests.fault_evaluation.golden import _coordinates, assert_matches, capture
from tests.fault_evaluation.models import build_model, sample_points

MODELS = [
    "f1_u0",
    "f4_u0",
    "f4_u3",
    "cutting_faults",
    "plugin",
    "stopped_fault",
    "faulted_fold",
    "onlap_domain",
]
ENV_VAR = "LOOPSTRUCTURAL_NATIVE_FAULTS"
# Map a model name to the set of feature names that cannot use the native path by design
# (a fold, or an interpolator that is not a DiscreteInterpolator on a StructuredGrid).
# Add an entry only with a comment that gives the reason.
FIELD_EXCLUSIONS: dict = {
    # The fold interpolator uses a tetrahedral mesh. Spec 007 supports only structured grids
    # (spec.md "Out of scope").
    "faulted_fold": {"folded_unit"},
}
_CACHE: dict = {}


def _case(name):
    """Return (model, points, box_size) and build each model once."""
    if name not in _CACHE:
        model = build_model(name)
        box = model.bounding_box
        size = float(np.max(np.asarray(box.maximum) - np.asarray(box.origin)))
        _CACHE[name] = (model, sample_points(name, 2000, seed=1), size)
    return _CACHE[name]


def _faulted(model):
    return [f for f in _coordinates(model) if getattr(f, "faults", None)]


def test_field_for_exists(monkeypatch):
    """Story 3: the native field of a faulted feature is available."""
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", True)
    features = _faulted(_case("f1_u0")[0])
    assert features
    for feature in features:
        assert _native_faults.field_for(feature) is not None


@pytest.mark.parametrize("name", MODELS)
def test_native_capture_matches_python(name, monkeypatch, record_property):
    """Values, gradients and model ids with native on equal those with native off."""
    model, points, size = _case(name)
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    # The native path must run, or this test compares the Python path with itself.
    assert _native_faults.native_enabled()
    for feature in _faulted(model):
        assert _native_faults.chain_for(feature) is not None
        if feature.name in FIELD_EXCLUSIONS.get(name, set()):
            continue
        assert _native_faults.field_for(feature) is not None
    native = capture(model, points)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
    reference = capture(model, points)
    # assert_matches lets model ids differ only where the Python value is within
    # 1e-8 * box_size of a unit boundary (FR-005b), using unit_boundary_distances.
    assert_matches(name, native, reference, size)
    worst = 0.0
    for key, value in reference.items():
        a, b = np.asarray(native[key]), np.asarray(value)
        if a.dtype.kind == "f" and a.shape == b.shape and np.isfinite(a).any():
            both = np.isfinite(a) & np.isfinite(b)
            if both.any():
                worst = max(worst, float(np.max(np.abs(a[both] - b[both]))))
    record_property("max_native_diff", worst)


def _point_on_fault(feature):
    """Find a point where the first fault frame value is 0, by bisection."""
    points = _case("f1_u0")[1]
    frame = feature.faults[0][0]
    side = frame.evaluate_value(points)
    a, b = points[np.argmax(side)], points[np.argmin(side)]
    assert frame.evaluate_value(a[None])[0] > 0 > frame.evaluate_value(b[None])[0]
    for _ in range(80):
        mid = 0.5 * (a + b)
        if frame.evaluate_value(mid[None])[0] > 0:
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


def test_gradient_on_fault_is_nan_and_warns_once(monkeypatch, caplog):
    """Story 3, scenario 2: a point on a fault has a NaN gradient and one warning."""
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    feature = _faulted(_case("f1_u0")[0])[0]
    point = _point_on_fault(feature)
    with caplog.at_level(logging.WARNING):
        gradient = feature.evaluate_gradient(point[None])
    assert gradient.shape == (1, 3)
    assert np.isnan(gradient).all()
    warnings = [r for r in caplog.records if "too close to a fault" in r.getMessage()]
    assert len(warnings) == 1
    # The Python path must also give a NaN gradient at this point.
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
    python_gradient = feature.evaluate_gradient(point[None])
    assert python_gradient.shape == (1, 3)
    assert np.isnan(python_gradient).all()


def test_unconformity_region_uses_native_path(monkeypatch):
    """A fault frame with an unconformity region still uses the native path."""
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    features = [
        f
        for f in _faulted(_case("onlap_domain")[0])
        if any(type(r).__name__ == "UnconformityFeature" for x in f.faults for r in x[0].regions)
    ]
    assert features
    for feature in features:
        assert _native_faults.chain_for(feature) is not None
