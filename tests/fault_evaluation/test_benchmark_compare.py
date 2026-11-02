"""Contract tests for benchmark result comparison (schema v1)."""

from __future__ import annotations

import copy
import importlib

import pytest


def _api(module_name: str, name: str):
    """Import a focused API and report an actionable absent-tool failure."""
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as error:
        if error.name == module_name or module_name.startswith(f"{error.name}."):
            pytest.fail(f"{module_name}.{name} is required by D2; benchmark compare is missing")
        raise
    value = getattr(module, name, None)
    assert value is not None, f"{module_name}.{name} is required by D2"
    return value


def _provenance(**overrides) -> dict:
    """Return a complete schema v1 provenance block with optional overrides."""
    provenance = {
        "python": "3.11.9",
        "numpy": "2.1.0",
        "scipy": "1.14.0",
        "loopstructural": "2.0.0",
        "os": "test",
        "cpu": "test cpu",
        "ram_bytes": 1_000_000,
    }
    provenance.update(overrides)
    return provenance


def _case(
    case_id: str,
    *,
    median: float = 2.0,
    state: str = "cold",
    points: int = 5000,
    model: str = "f1_u0",
    operation: str = "oldest_value",
) -> dict:
    """Return one complete schema v1 case row with a fixed median time."""
    return {
        "case_id": case_id,
        "model": model,
        "operation": operation,
        "points": points,
        "state": state,
        "rounds": 2,
        "time_samples_seconds": [median, median],
        "median_seconds": median,
        "fault_calls_total": 2,
        "fault_calls_top_level": 1,
        "field_calls_total": 3,
        "fault_points_total": 4,
        "field_points_total": 7,
        "peak_extra_rss_bytes": 4096,
        "rss_unit": "bytes",
        "code_revision": "abc123",
        "dirty": False,
        "prerequisite_parent_revision": "parent",
        "prerequisite_fix_revision": "fix",
        "baseline_revision": "baseline",
        "input_hash": "0" * 64,
    }


def _result(cases: list[dict], **provenance_overrides) -> dict:
    """Return a complete schema v1 result document for the given rows."""
    return {
        "schema_version": 1,
        "pass": "merged",
        "cases": cases,
        "provenance": _provenance(**provenance_overrides),
        "pytest_benchmark": {"raw": {"stats": {"mean": 2.0}}},
    }


def _row(result: dict, case_id: str) -> dict:
    """Return the comparison row for one case id."""
    rows = {row["case_id"]: row for row in result["cases"]}
    return rows[case_id]


def test_ratio_and_slower_mark_for_a_two_to_three_second_change():
    """A 2 s to 3 s change gives ratio 1.5 and is marked slower."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold", median=2.0)])
    new = _result([_case("tiny:cold", median=3.0)])
    result = compare_results(base, new)
    row = _row(result, "tiny:cold")
    assert row["ratio"] == pytest.approx(1.5)
    assert row["slower"] is True


def test_small_change_within_ten_percent_is_not_marked_slower():
    """A 2 s to 2.1 s change (ratio 1.05) is within the 10 percent noise band."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold", median=2.0)])
    new = _result([_case("tiny:cold", median=2.1)])
    result = compare_results(base, new)
    row = _row(result, "tiny:cold")
    assert row["ratio"] == pytest.approx(1.05)
    assert row["slower"] is False


def test_zero_baseline_gives_explicit_undefined_ratio():
    """A zero baseline time gives no numeric ratio and an explicit undefined status."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold", median=0.0)])
    new = _result([_case("tiny:cold", median=1.0)])
    result = compare_results(base, new)
    row = _row(result, "tiny:cold")
    assert row["ratio"] is None
    assert row["ratio_status"] == "undefined_zero_baseline"


@pytest.mark.parametrize(
    "base_cases,new_cases,match",
    [
        pytest.param(
            [_case("tiny:cold"), _case("other:cold")],
            [_case("tiny:cold")],
            "missing",
            id="case-missing-from-new",
        ),
        pytest.param(
            [_case("tiny:cold")],
            [_case("tiny:cold"), _case("tiny:cold")],
            "duplicate",
            id="duplicate-case-in-new",
        ),
        pytest.param(
            [_case("tiny:cold", points=5000)],
            [_case("tiny:cold", points=50000)],
            "points",
            id="different-point-count",
        ),
        pytest.param(
            [_case("tiny:cold", state="cold")],
            [_case("tiny:cold", state="warm")],
            "state",
            id="cold-baseline-warm-new",
        ),
    ],
)
def test_incompatible_or_missing_cases_raise_value_error(base_cases, new_cases, match):
    """Missing, duplicate, resized or cold/warm-mixed cases are rejected, never skipped."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result(base_cases)
    new = _result(new_cases)
    with pytest.raises(ValueError, match=match):
        compare_results(base, new)


def test_duplicate_case_in_base_raises_value_error():
    """A duplicate case id in the baseline document is rejected."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold"), _case("tiny:cold")])
    new = _result([_case("tiny:cold")])
    with pytest.raises(ValueError, match="duplicate"):
        compare_results(base, new)


def test_fault_and_unconformity_scaling_tables_are_built_from_model_names():
    """Scaling tables group ratios by fault count (F) and unit count (U)."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    models = ["f0_u0", "f1_u0", "f2_u0", "f0_u1"]
    base = _result([_case(f"{model}:cold", model=model) for model in models])
    new = _result([_case(f"{model}:cold", model=model, median=3.0) for model in models])
    result = compare_results(base, new)
    scaling = result["scaling"]
    assert set(scaling) >= {"fault", "unconformity"}
    # Fault table has one row per fault count (0, 1, 2) at unit count 0.
    assert len(scaling["fault"]) == 3
    # Unconformity table has one row per unit count (0, 1) at fault count 0.
    assert len(scaling["unconformity"]) == 2
    assert all(row["ratio"] == pytest.approx(1.5) for row in scaling["fault"])


def test_provenance_difference_produces_a_warning():
    """A changed Python version is reported as a provenance warning."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold")])
    new = _result([_case("tiny:cold")], python="3.12.4")
    result = compare_results(base, new)
    assert any(
        "provenance" in text.lower() and "python" in text.lower() for text in result["warnings"]
    )


def test_identical_provenance_produces_no_provenance_warning():
    """Matching provenance blocks produce no provenance warning."""
    compare_results = _api("benchmarks.fault_evaluation.compare", "compare_results")
    base = _result([_case("tiny:cold")])
    new = copy.deepcopy(base)
    result = compare_results(base, new)
    assert not [text for text in result["warnings"] if "provenance" in text.lower()]
