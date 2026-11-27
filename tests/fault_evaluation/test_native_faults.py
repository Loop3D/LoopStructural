"""Native restoration agrees with the Python path (spec 007, story 2, FR-003, FR-005).

The native result comes from ``_native_faults.restore``. The reference comes from
``restore_through`` with the native path turned off.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("loop_interpolation_core")

from LoopStructural import LoopStructuralConfig
from LoopStructural.modelling.features import _native_faults
from LoopStructural.utils._fault_block import check_block, restore_through
from tests.fault_evaluation.golden import _coordinates
from tests.fault_evaluation.models import build_model, sample_points

MODELS = ["f1_u0", "f4_u0", "f4_u3", "cutting_faults", "plugin", "stopped_fault", "faulted_fold"]
MUST_BE_SUPPORTED = ["f1_u0", "f4_u0", "f4_u3", "cutting_faults"]
ENV_VAR = "LOOPSTRUCTURAL_NATIVE_FAULTS"

_CACHE: dict = {}


def _case(name):
    """Return (model, points, box_size, faulted features) and build each model once."""
    if name not in _CACHE:
        model = build_model(name)
        points = sample_points(name, 2000, seed=1)
        # Same rule as ``_box_size`` in regenerate.py: the largest box extent.
        box = model.bounding_box
        size = float(np.max(np.asarray(box.maximum) - np.asarray(box.origin)))
        features = [
            f
            for f in _coordinates(model)
            if getattr(f, "faults", None) and getattr(f, "faults_enabled", False)
        ]
        _CACHE[name] = (model, points, size, features)
    return _CACHE[name]


@pytest.fixture
def python_path(monkeypatch):
    """Turn the native path off for the reference and clear the variable."""
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)


def _variants(feature):
    """Return the keyword arguments of each call to compare."""
    first = feature.faults[0]
    cases = [{}, {"reverse": True}]
    for fault in feature.faults:
        cases += [{"block": {fault.name: side}} for side in (1, -1, 0)]
    if len(feature.faults) >= 2:
        cases.append({"restored_through": (first,)})
    return cases


@pytest.mark.parametrize("name", MODELS)
def test_native_restore_matches_python(name, monkeypatch, record_property):
    """Native and Python restoration agree for every supported faulted feature."""
    _, points, size, features = _case(name)
    skipped, supported, worst = [], [], 0.0
    for feature in features:
        monkeypatch.delenv(ENV_VAR, raising=False)
        monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
        if _native_faults.chain_for(feature) is None:
            skipped.append(feature.name)
            continue
        supported.append(feature.name)
        for kwargs in _variants(feature):
            monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
            native = _native_faults.restore(feature, points, **kwargs)
            assert native is not None
            monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
            reference = restore_through(
                points,
                feature.faults,
                reverse=kwargs.get("reverse", False),
                block=check_block(kwargs.get("block")),
                restored_through=kwargs.get("restored_through", ()),
                owner=None,
            )
            assert native.shape == reference.shape
            # NaN rows must be the same rows (FR-005a).
            nan = np.isnan(native).any(axis=1)
            assert np.array_equal(nan, np.isnan(reference).any(axis=1)), (feature.name, kwargs)
            if not nan.all():
                diff = float(np.nanmax(np.abs(native - reference)))
                worst = max(worst, diff)
                assert diff <= 1e-8 * size, (feature.name, kwargs, diff)
    record_property("max_native_diff", worst)
    record_property("native_supported_features", ",".join(supported))
    record_property("native_skipped_features", ",".join(skipped))
    if name in MUST_BE_SUPPORTED:
        assert features
        assert not skipped


def _supported_feature():
    return _case("f1_u0")[3][0]


def test_native_used_by_default(monkeypatch):
    """Story 4, scenario 1: with the defaults the native path is used."""
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    points = _case("f1_u0")[1]
    assert isinstance(_native_faults.restore(_supported_feature(), points), np.ndarray)


def test_native_off_by_setting(monkeypatch):
    """Story 4, scenario 2: the setting False turns the native path off."""
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
    points = _case("f1_u0")[1]
    assert _native_faults.restore(_supported_feature(), points) is None


def test_native_off_by_variable(monkeypatch):
    """Story 4, scenario 3: the variable "0" turns it off; the setting wins over it."""
    points = _case("f1_u0")[1]
    feature = _supported_feature()
    monkeypatch.setenv(ENV_VAR, "0")
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    assert _native_faults.restore(feature, points) is None
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", True)
    assert isinstance(_native_faults.restore(feature, points), np.ndarray)
