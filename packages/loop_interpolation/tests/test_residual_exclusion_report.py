import dataclasses

import numpy as np
import pytest
from loop_interpolation._diagnostics import (
    UNKNOWN_EXPLANATION,
    ConstraintDiagnosticsReport,
    ConstraintFamilyDiagnostics,
    ExclusionRecord,
    RegionCoverageDiagnostics,
    ResidualSummary,
    build_exclusion_records,
    direction_cosines,
    summarise_residuals,
)


def _family(
    name,
    dropped_rows=0,
    source_point_count=0,
    outside_model_point_count=0,
    row_count=10,
):
    return ConstraintFamilyDiagnostics(
        name=name,
        active=True,
        row_count=row_count,
        dropped_rows=dropped_rows,
        effective_weight_mean=1.0,
        effective_weight_min=1.0,
        effective_weight_max=1.0,
        source_point_count=source_point_count,
        outside_model_point_count=outside_model_point_count,
    )


def _unattributed(records):
    return [r for r in records if r.reason == "unattributed"]


def test_direction_cosines_pure_function():
    model = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    data = np.array(
        [
            [2.0, 0.0, 0.0],
            [-3.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [1.0, 0.0, 0.0],
            [np.nan, 0.0, 0.0],
        ]
    )
    model_before = model.copy()
    data_before = data.copy()
    result = direction_cosines(model, data)
    assert result.shape == (5,)
    np.testing.assert_allclose(result[:3], [1.0, -1.0, 0.0], atol=1e-12)
    assert np.isnan(result[3])
    assert np.isnan(result[4])
    # Inputs are not changed.
    np.testing.assert_array_equal(model, model_before)
    np.testing.assert_array_equal(data, data_before)
    # A shape mismatch is an error.
    with pytest.raises(ValueError):
        direction_cosines(model[:3], data[:2])
    with pytest.raises(ValueError):
        direction_cosines(model, data[:, :2])
    # The result does not depend on the length of the vectors.
    a = np.array([[1.0, 2.0, 3.0], [0.0, 1.0, 1.0]])
    b = np.array([[3.0, 1.0, -2.0], [1.0, 1.0, 0.0]])
    np.testing.assert_allclose(
        direction_cosines(a, b), direction_cosines(5.0 * a, 0.25 * b), atol=1e-12
    )


def test_report_without_verification_is_unchanged():
    family = _family("value", source_point_count=4)
    region = RegionCoverageDiagnostics(
        total_support_nodes=10,
        active_region_nodes=6,
        inactive_region_nodes=4,
        active_fraction=0.6,
    )
    report = ConstraintDiagnosticsReport(
        interpolator_type="FiniteDifferenceInterpolator",
        families={"value": family},
        region_coverage=region,
        outside_model_points={"value": 1},
    )
    data = report.to_dict()
    assert set(data.keys()) == {
        "interpolator_type",
        "total_rows",
        "families",
        "region_coverage",
        "outside_model_points",
    }
    assert data["total_rows"] == 10  # row_count of the one family
    assert set(data["families"]["value"].keys()) == {
        "active",
        "row_count",
        "dropped_rows",
        "effective_weight_mean",
        "effective_weight_min",
        "effective_weight_max",
        "source_point_count",
        "outside_model_point_count",
    }
    assert set(data["region_coverage"].keys()) == {
        "total_support_nodes",
        "active_region_nodes",
        "inactive_region_nodes",
        "active_fraction",
    }
    # The old fields keep their order. New fields can only come after them.
    assert [f.name for f in dataclasses.fields(ConstraintDiagnosticsReport)][:4] == [
        "interpolator_type",
        "families",
        "region_coverage",
        "outside_model_points",
    ]
    assert [f.name for f in dataclasses.fields(ConstraintFamilyDiagnostics)][:9] == [
        "name",
        "active",
        "row_count",
        "dropped_rows",
        "effective_weight_mean",
        "effective_weight_min",
        "effective_weight_max",
        "source_point_count",
        "outside_model_point_count",
    ]


def test_builder_records_keep_their_own_key():
    families = {"restricted_orthogonal": _family("restricted_orthogonal", dropped_rows=3)}
    before = dataclasses.replace(families["restricted_orthogonal"])
    records = build_exclusion_records(
        {"restricted_orthogonal": {"invalid_reference_direction": 3}}, families
    )
    assert len(records) == 1
    record = records[0]
    assert isinstance(record, ExclusionRecord)
    assert record.family == "restricted_orthogonal"
    assert record.reason == "invalid_reference_direction"
    assert record.count == 3
    assert _unattributed(records) == []
    assert families["restricted_orthogonal"] == before


def test_unknown_reason_keeps_its_count():
    families = {"value": _family("value", dropped_rows=2)}
    records = build_exclusion_records({"value": {"x": 2}}, families)
    assert len(records) == 1
    assert records[0].reason == "x"
    assert records[0].count == 2
    assert records[0].explanation == UNKNOWN_EXPLANATION


def test_unattributed_rows_are_reported():
    families = {
        "value": _family(
            "value", dropped_rows=5, source_point_count=10, outside_model_point_count=2
        )
    }
    records = build_exclusion_records({}, families)
    unattributed = _unattributed(records)
    assert len(unattributed) == 1
    assert unattributed[0].family == "value"
    assert unattributed[0].count == 3


def test_gradient_outside_points_give_no_false_unattributed():
    def make(dropped):
        return {
            "gradient": _family(
                "gradient",
                dropped_rows=dropped,
                source_point_count=10,
                outside_model_point_count=2,
            )
        }

    records = build_exclusion_records({}, make(4), rows_per_point={"gradient": 2})
    assert records == ()
    records = build_exclusion_records({}, make(6), rows_per_point={"gradient": 2})
    unattributed = _unattributed(records)
    assert len(unattributed) == 1
    assert unattributed[0].count == 2


def test_dropped_rows_none_gives_no_unattributed():
    families = {"interface": _family("interface", dropped_rows=None)}
    assert build_exclusion_records({}, families) == ()


def test_surfe_style_dropped_equals_outside():
    families = {
        "value": _family(
            "value", dropped_rows=2, source_point_count=10, outside_model_point_count=2
        )
    }
    assert build_exclusion_records({}, families, rows_per_point=None) == ()


def test_residual_summary_is_importable():
    # The name must exist. Its behaviour is tested with the implementation task.
    assert ResidualSummary is not None
    assert callable(summarise_residuals)
