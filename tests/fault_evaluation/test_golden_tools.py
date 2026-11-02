"""Contract tests for golden capture, comparison, and explicit regeneration."""

from __future__ import annotations

import importlib
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


def _api(module_name: str, function_name: str):
    """Import a required golden helper and give a focused missing-API failure."""
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as error:
        if error.name != module_name:
            raise
        pytest.fail(f"{module_name}.{function_name} API is required by D1")
    function = getattr(module, function_name, None)
    assert callable(function), f"{module_name}.{function_name} API is required by D1"
    return function


def _pair(**updates):
    """Return a small captured-output mapping with one point and one unit boundary."""
    actual = {
        "model_ids": np.array([2], dtype=np.int64),
        "feature::frame::value": np.array([3.0]),
        "feature::frame::gradient": np.array([[1.0, 2.0, 3.0]]),
        "restored::fault_0>fault_1": np.array([[1.0, 2.0, 3.0]]),
        "surface::unit@50::vertex_count": np.array([4]),
        "surface::unit@50::bounds": np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
        "surface::unit@50::area": np.array([2.0]),
        # Distances are derived from declared scalar unit bounds. Each row is a
        # sample and each column is one boundary; +inf marks an ineligible value.
        "unit_boundary_distances": np.array([[0.0, 4.0]]),
    }
    actual.update(updates)
    expected = {key: value.copy() for key, value in actual.items()}
    return actual, expected


def test_capture_has_flat_numeric_keys_for_all_d1_output_families(monkeypatch):
    capture = _api("tests.fault_evaluation.golden", "capture")
    # This solved special case exercises parent and frame-coordinate evaluation
    # without invoking a solver. Replace surface extraction with one tiny mesh.
    from tests.fault_evaluation.models import build_model, sample_points

    points = sample_points("cutting_faults", 2, seed=17)
    model = build_model("cutting_faults")
    mesh = SimpleNamespace(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        triangles=np.array([[0, 1, 2]], dtype=np.int64),
    )
    for feature in model.features:
        monkeypatch.setattr(
            feature,
            "surfaces",
            lambda values, bounding_box, **kwargs: [mesh for _ in values],
        )
    captured = capture(model, points)

    assert isinstance(captured, dict)
    assert captured
    assert all(isinstance(key, str) for key in captured)
    assert all(isinstance(value, np.ndarray) for value in captured.values())
    assert all(np.issubdtype(value.dtype, np.number) for value in captured.values())
    assert "model_ids" in captured
    assert "unit_boundary_distances" in captured
    assert captured["model_ids"].shape == (len(points),)
    assert captured["unit_boundary_distances"].ndim == 2
    assert captured["unit_boundary_distances"].shape[0] == len(points)
    assert captured["unit_boundary_distances"].shape[1] > 0
    coordinates = {
        coordinate.name
        for feature in model.features
        for coordinate in (list(getattr(feature, "features", None) or []) + [feature])
        if getattr(coordinate, "builder", None) is not None
    }
    assert coordinates
    for name in coordinates:
        assert f"feature::{name}::value" in captured
        assert f"feature::{name}::gradient" in captured
    assert "restored::fault_1>fault_0" in captured
    assert any(key.startswith("surface::") and key.endswith("::vertex_count") for key in captured)
    assert any(key.startswith("surface::") and key.endswith("::bounds") for key in captured)
    assert any(key.startswith("surface::") and key.endswith("::area") for key in captured)
    summary = next(key.rsplit("::", 1)[0] for key in captured if key.endswith("::vertex_count"))
    np.testing.assert_array_equal(captured[f"{summary}::vertex_count"], [3])
    np.testing.assert_allclose(captured[f"{summary}::bounds"], [[0, 0, 0], [1, 1, 0]])
    np.testing.assert_allclose(captured[f"{summary}::area"], [0.5])


def test_capture_records_empty_surfaces_when_all_levels_are_outside_range(monkeypatch):
    from tests.fault_evaluation.golden import capture
    from tests.fault_evaluation.models import build_model, sample_points

    model = build_model("f0_u0")
    group = model.stratigraphic_column.get_groups()[0]
    for index, unit in enumerate(group.units):
        unit.min_value, unit.max_value = 2000 + index * 100, 2100 + index * 100
    feature = next(feature for feature in model.features if feature.name == group.name)
    monkeypatch.setattr(feature, "surfaces", lambda values, bounding_box, **kwargs: [])
    captured = capture(model, sample_points("f0_u0", 2))
    summaries = [key.rsplit("::", 1)[0] for key in captured if key.endswith("::vertex_count")]
    assert summaries
    for summary in summaries:
        np.testing.assert_array_equal(captured[f"{summary}::vertex_count"], [0])
        np.testing.assert_array_equal(captured[f"{summary}::bounds"], np.zeros((2, 3)))
        np.testing.assert_array_equal(captured[f"{summary}::area"], [0.0])


def test_capture_rejects_unexpected_empty_surfaces(monkeypatch):
    from tests.fault_evaluation.golden import capture
    from tests.fault_evaluation.models import build_model, sample_points

    model = build_model("f0_u0")
    group = model.stratigraphic_column.get_groups()[0]
    feature = next(feature for feature in model.features if feature.name == group.name)
    monkeypatch.setattr(feature, "surfaces", lambda values, bounding_box, **kwargs: [])
    with pytest.raises(ValueError, match="Empty surface"):
        capture(model, sample_points("f0_u0", 2))


def test_comparator_accepts_box_scaled_tolerance_and_boundary_id_only():
    compare = _api("tests.fault_evaluation.golden", "assert_matches")
    actual, expected = _pair()
    actual["feature::frame::value"][0] += 1e-7
    actual["model_ids"][0] = 7
    # The value perturbation is within 1e-8 * box_size, and this point is
    # exactly on a declared scalar boundary.
    compare("tiny", actual, expected, box_size=100.0)

    actual, expected = _pair()
    actual["model_ids"][0] = 7
    actual["unit_boundary_distances"][:] = 1e-3
    expected["unit_boundary_distances"][:] = 1e-3
    with pytest.raises(AssertionError, match=r"tiny.*model_ids.*max.*point"):
        compare("tiny", actual, expected, box_size=100.0)

    actual, expected = _pair()
    actual["model_ids"][0] = 7
    actual["unit_boundary_distances"][:] = np.inf
    expected["unit_boundary_distances"][:] = np.inf
    with pytest.raises(AssertionError, match="model_ids"):
        compare("tiny", actual, expected, box_size=100.0)


@pytest.mark.parametrize(
    ("key", "replacement"),
    [
        ("feature::frame::value", np.array([3.01])),
        ("feature::frame::gradient", np.array([[1.01, 2.0, 3.0]])),
        ("restored::fault_0>fault_1", np.array([[1.0, 2.01, 3.0]])),
        ("surface::unit@50::bounds", np.array([[0.0, 0.0, 0.0], [1.01, 1.0, 1.0]])),
        ("surface::unit@50::area", np.array([2.01])),
    ],
)
def test_comparator_reports_model_output_max_difference_and_point(key, replacement):
    compare = _api("tests.fault_evaluation.golden", "assert_matches")
    actual, expected = _pair()
    actual[key] = replacement
    with pytest.raises(AssertionError) as error:
        compare("f0_u0", actual, expected, box_size=1.0)
    message = str(error.value)
    assert "f0_u0" in message
    assert key in message
    assert "max" in message.lower()
    assert "point" in message.lower()


def test_comparator_rejects_mismatched_keys_shapes_and_nonfinite_masks():
    compare = _api("tests.fault_evaluation.golden", "assert_matches")
    actual, expected = _pair()
    actual["extra"] = np.array([1.0])
    with pytest.raises(AssertionError, match="key"):
        compare("case", actual, expected, box_size=1.0)

    actual, expected = _pair()
    actual["feature::frame::gradient"] = np.ones((1, 2))
    with pytest.raises(AssertionError, match="shape"):
        compare("case", actual, expected, box_size=1.0)

    actual, expected = _pair()
    actual["feature::frame::value"][0] = np.nan
    with pytest.raises(AssertionError, match="NaN|finite|mask"):
        compare("case", actual, expected, box_size=1.0)


def test_comparator_requires_exact_surface_counts_and_zero_stays_zero():
    compare = _api("tests.fault_evaluation.golden", "assert_matches")
    actual, expected = _pair()
    actual["surface::unit@50::vertex_count"][0] += 1
    with pytest.raises(AssertionError, match="vertex_count"):
        compare("surface-case", actual, expected, box_size=1.0)

    actual, expected = _pair()
    actual["surface::unit@50::area"][0] = 1e-12
    expected["surface::unit@50::area"][0] = 0.0
    with pytest.raises(AssertionError, match="area"):
        compare("surface-case", actual, expected, box_size=1.0)


def test_comparator_rejects_mismatching_fault_restoration():
    compare = _api("tests.fault_evaluation.golden", "assert_matches")
    actual, expected = _pair()
    actual["restored::fault_0>fault_1"][0, 0] += 0.01
    with pytest.raises(AssertionError, match="restored::fault_0>fault_1"):
        compare("fault-case", actual, expected, box_size=1.0)


def test_missing_reference_and_validation_mode_do_not_write(tmp_path, monkeypatch):
    regenerate_main = _api("tests.fault_evaluation.regenerate", "main")
    regenerate = importlib.import_module("tests.fault_evaluation.regenerate")
    monkeypatch.setattr(regenerate, "REFERENCE_DIR", tmp_path, raising=False)
    monkeypatch.setattr(regenerate, "model_names", lambda: ["fake"], raising=False)
    monkeypatch.setattr(regenerate, "build_model", lambda name: object(), raising=False)
    monkeypatch.setattr(
        regenerate, "sample_points", lambda *args, **kwargs: np.zeros((1, 3)), raising=False
    )
    monkeypatch.setattr(
        regenerate, "capture", lambda model, points: {"model_ids": np.array([0])}, raising=False
    )

    before = sorted(path.name for path in tmp_path.iterdir())
    result = regenerate_main([])
    after = sorted(path.name for path in tmp_path.iterdir())
    assert result == 1
    assert before == after
    assert not (tmp_path / "fake.npz").exists()


def test_write_requires_explicit_flag_and_targets_only_requested_directory(tmp_path, monkeypatch):
    regenerate_main = _api("tests.fault_evaluation.regenerate", "main")
    regenerate = importlib.import_module("tests.fault_evaluation.regenerate")
    output = tmp_path / "references"
    output.mkdir()
    fake_model = SimpleNamespace(
        features=[], bounding_box=SimpleNamespace(origin=np.zeros(3), maximum=np.ones(3))
    )
    monkeypatch.setattr(regenerate, "REFERENCE_DIR", output, raising=False)
    monkeypatch.setattr(regenerate, "model_names", lambda: ["fake"], raising=False)
    monkeypatch.setattr(regenerate, "build_model", lambda name: fake_model, raising=False)
    monkeypatch.setattr(
        regenerate, "sample_points", lambda *args, **kwargs: np.zeros((1, 3)), raising=False
    )
    monkeypatch.setattr(
        regenerate, "capture", lambda model, points: {"model_ids": np.array([0])}, raising=False
    )

    result = regenerate_main(["--write"])
    assert result == 0
    assert (output / "fake.npz").exists()
    assert (output / "manifest.json").exists()
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf8"))
    archive_bytes = (output / "fake.npz").read_bytes()
    assert manifest["models"]["fake"]["sha256"] == hashlib.sha256(archive_bytes).hexdigest()
    for revision in (
        "baseline_revision",
        "prerequisite_fix_revision",
        "prerequisite_parent_revision",
    ):
        assert isinstance(manifest["provenance"][revision], str)
        assert manifest["provenance"][revision]
    with np.load(output / "fake.npz") as archive:
        np.testing.assert_array_equal(archive["model_ids"], [0])
    changed = {path.relative_to(tmp_path) for path in tmp_path.rglob("*") if path.is_file()}
    assert changed == {Path("references/fake.npz"), Path("references/manifest.json")}
