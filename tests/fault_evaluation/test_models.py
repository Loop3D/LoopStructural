"""Check the fixed model catalogue and its stored fields."""

from __future__ import annotations

import numpy as np
import pytest

SYNTHETIC_NAMES = [
    f"f{faults}_u{unconformities}" for faults in (0, 1, 2, 4) for unconformities in range(4)
]
SPECIAL_NAMES = ["stopped_fault", "cutting_faults", "onlap_domain", "faulted_fold", "plugin"]
EXPECTED_NAMES = SYNTHETIC_NAMES + SPECIAL_NAMES


def _api(name):
    """Load a model helper and report its public test contract when it is absent."""
    import importlib

    try:
        models = importlib.import_module("tests.fault_evaluation.models")
    except ModuleNotFoundError as error:
        if error.name != "tests.fault_evaluation.models":
            raise
        pytest.fail(
            f"tests.fault_evaluation.models.{name} API is required by D1 (models module is missing)"
        )

    helper = getattr(models, name, None)
    assert callable(helper), f"tests.fault_evaluation.models.{name} API is required by D1"
    return helper


def _builders(feature):
    builder = feature.builder
    return list(getattr(builder, "builders", None) or [builder])


def _feature_list(model):
    return list(model.features)


def test_model_names_matches_d1_catalogue():
    model_names = _api("model_names")
    assert model_names() == EXPECTED_NAMES, "model_names() must return the 21 D1 names in order"


@pytest.mark.parametrize("name", EXPECTED_NAMES)
def test_build_model_relationships_and_stored_fields(name, monkeypatch):
    """Build each fixture without invoking any solver."""
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator
    from loop_interpolation._geological_interpolator import GeologicalInterpolator

    calls = []

    def reject_solve(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError(f"build_model({name!r}) must load solved fields without solver calls")

    monkeypatch.setattr(DiscreteInterpolator, "solve_system", reject_solve)
    monkeypatch.setattr(GeologicalInterpolator, "solve_system", reject_solve)
    build_model = _api("build_model")
    model = build_model(name)
    assert not calls, f"build_model({name!r}) called a solver"

    features = _feature_list(model)
    assert features, f"build_model({name!r}) returned no model features"
    assert model.stratigraphic_column is not None, f"{name} must include a stratigraphic column"
    assert (
        len(model.stratigraphic_column.order) > 0
    ), f"{name} stratigraphic column must not be empty"

    for feature in features:
        for builder in _builders(feature):
            assert builder._up_to_date, f"{name}/{feature.name} stored builder state is not current"
            assert (
                builder.interpolator.up_to_date
            ), f"{name}/{feature.name} stored field is not current"

    if name in SYNTHETIC_NAMES:
        fault_count = int(name.split("_u")[0][1:])
        unconformity_count = int(name.rsplit("_u", 1)[1])
        faults = [feature for feature in features if getattr(feature, "type", None).name == "FAULT"]
        groups = [model.get_feature_by_name(f"group_{i}") for i in range(unconformity_count + 1)]
        assert len(faults) == fault_count, f"{name} fault catalogue is incomplete"
        assert len(groups) == unconformity_count + 1, f"{name} group catalogue is incomplete"
        # The data recipe creates group_0 first (youngest) and the last group is oldest.
        group_positions = [features.index(group) for group in groups]
        assert group_positions == sorted(
            group_positions
        ), f"{name} group creation age order is incorrect"
        expected_fault_order = [f"fault_{i}" for i in reversed(range(fault_count))]
        for group in groups:
            assert [
                fault.name for fault in group.faults
            ] == expected_fault_order, (
                f"{name}/{group.name} fault order differs from the baseline recipe"
            )
        for index in range(unconformity_count):
            younger = groups[index]
            assert any(
                getattr(region, "parent", None) is younger for region in groups[index + 1].regions
            ), f"{name}/{groups[index + 1].name} must be bounded by the preceding younger group"

    if name == "stopped_fault":
        group_fault_sets = {
            tuple(f.name for f in feature.faults)
            for feature in features
            if feature.name.startswith("group_")
        }
        assert (
            len(group_fault_sets) >= 2
        ), "stopped_fault must have different fault sets across its boundary"
    elif name == "cutting_faults":
        faults = [feature for feature in features if getattr(feature, "type", None).name == "FAULT"]
        assert len(faults) >= 2, "cutting_faults must contain at least two faults"
        assert any(
            any(parent is other for parent in fault.faults) for fault in faults for other in faults
        ), "cutting_faults must preserve the fault-on-fault relationship"
    elif name == "onlap_domain":
        from LoopStructural.modelling.features import FeatureType

        assert any(
            region.type == FeatureType.ONLAPUNCONFORMITY
            for feature in features
            for region in feature.regions
        ), "onlap_domain must contain an onlap unconformity region"
        assert any(
            getattr(feature, "type", None).name == "DOMAINFAULT"
            or getattr(feature, "is_domain_fault", False)
            for feature in features
        ), "onlap_domain must contain a domain fault"
    elif name == "faulted_fold":
        from LoopStructural.modelling.features.fold import FoldFrame

        folded = [feature for feature in features if getattr(feature, "fold", None) is not None]
        frames = [feature for feature in features if isinstance(feature, FoldFrame)]
        assert folded, "faulted_fold must contain a feature linked to a FoldEvent"
        assert frames, "faulted_fold must contain a FoldFrame"
        assert any(frame.faults for frame in frames), "faulted_fold frame must be cut by a fault"
    elif name == "plugin":
        from pathlib import Path

        data_dir = (
            Path(__file__).parent.parent.parent
            / "specs"
            / "003-fault-evaluation-performance"
            / "data"
            / "plugin_model"
        )
        reference = np.load(data_dir / "qgis_check.npz")
        points = reference["points"]
        group = model.get_feature_by_name("Group_0")
        actual = {"Group_0": group.evaluate_value(points)}
        for feature_name in ("7439", "2997", "3496", "3498", "12647", "12658"):
            actual[f"{feature_name}__0"] = model.get_feature_by_name(feature_name)[
                0
            ].evaluate_value(points)
        actual["restored_Group_0"] = group._apply_faults(points.copy())
        for key, value in actual.items():
            expected = reference[key]
            assert np.array_equal(
                np.isnan(value), np.isnan(expected)
            ), f"plugin {key} NaN mask differs from QGIS"
            np.testing.assert_allclose(
                value,
                expected,
                rtol=0.0,
                atol=1e-8,
                equal_nan=True,
                err_msg=f"plugin {key} differs from qgis_check.npz",
            )


@pytest.mark.parametrize("name", EXPECTED_NAMES)
def test_sample_points_is_fixed_and_inside_model_box(name):
    sample_points = _api("sample_points")
    first = sample_points(name, 64)
    second = sample_points(name, 64)
    assert isinstance(first, np.ndarray), f"sample_points({name!r}) must return an ndarray"
    assert first.shape == (64, 3), f"sample_points({name!r}) must return shape (count, 3)"
    assert np.isfinite(first).all(), f"sample_points({name!r}) must return finite points"
    np.testing.assert_array_equal(
        first, second, err_msg=f"sample_points({name!r}) must be deterministic"
    )
    model = _api("build_model")(name)
    box = model.bounding_box
    assert np.all(first >= box.origin) and np.all(
        first <= box.maximum
    ), f"sample_points({name!r}) must stay inside the model bounding box"
