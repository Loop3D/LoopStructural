"""Contract tests for numerical performance acceptance reporting (schema v1).

The acceptance check reads a base and a new schema v1 result document. It
reports one row for each target SC-002 to SC-010. A failed target is a
failure of the whole result. It is never a warning.
"""

from __future__ import annotations

import copy
import importlib

import pytest

ACCEPTANCE_MODULE = "benchmarks.fault_evaluation.acceptance"
GRADIENT_RULE = "op1-one-restoration-per-fault-set-and-point-batch"
CRITERIA = (
    "SC-002",
    "SC-003",
    "SC-004",
    "SC-005",
    "SC-006",
    "SC-007",
    "SC-008",
    "SC-009",
    "SC-010",
)
LIMITS = {
    "SC-002": 0.40,
    "SC-003": 0.15,
    "SC-004": 1.5,
    "SC-005": 0.50,
    "SC-006": 1.20,
    "SC-007": 1.10,
    "SC-008": 900.0,
    "SC-009": 0.50,
    "SC-010": 0.10,
}

# Case ids of the cases that the targets use. The id format matches run.py.
OLD_F4U3 = "f4_u3:oldest_value:5000:cold"
ALL_F4U3 = "f4_u3:all_group_values:5000:cold"
ALL_F4U0 = "f4_u0:all_group_values:5000:cold"
WARM_F4U3 = "f4_u3:all_group_values:5000:warm"
PLUGIN_OLD = "plugin:oldest_value:20000:cold"
SINGLE_F1U0 = "f1_u0:oldest_value:20000:cold"

# Row fields: model, operation, points, state, median seconds, fault calls,
# fault points, peak extra RSS bytes.
BASE_ROWS = [
    ("f4_u3", "oldest_value", 5000, "cold", 2.5, 16, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "cold", 8.65, 52, 5000, 1_000_000),
    ("f4_u0", "all_group_values", 5000, "cold", 0.67, 12, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "warm", 8.0, 52, 5000, 1_000_000),
    ("plugin", "oldest_value", 20000, "cold", 8.2, 6, 20000, 1_000_000),
    ("f1_u0", "oldest_value", 20000, "cold", 0.68, 2, 20000, 500_000),
]
NEW_ROWS = [
    ("f4_u3", "oldest_value", 5000, "cold", 1.0, 4, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "cold", 1.0, 4, 5000, 1_000_000),
    ("f4_u0", "all_group_values", 5000, "cold", 0.8, 4, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "warm", 0.05, 0, 0, 1_000_000),
    ("plugin", "oldest_value", 20000, "cold", 3.0, 6, 20000, 1_000_000),
    ("f1_u0", "oldest_value", 20000, "cold", 0.2, 2, 20000, 500_000),
]


def _api(module_name: str, name: str):
    """Import a focused API and report an actionable absent-tool failure."""
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as error:
        if error.name == module_name or module_name.startswith(f"{error.name}."):
            pytest.fail(f"{module_name}.{name} is required by D2; acceptance reporting is missing")
        raise
    value = getattr(module, name, None)
    assert value is not None, f"{module_name}.{name} is required by D2"
    return value


def _check_acceptance():
    """Return the check_acceptance API or fail with a message that names it."""
    return _api(ACCEPTANCE_MODULE, "check_acceptance")


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


def _case(row: tuple, code_revision: str) -> dict:
    """Return one complete schema v1 case row from a fixed test row."""
    model, operation, points, state, median, calls, fault_points, rss = row
    return {
        "case_id": f"{model}:{operation}:{points}:{state}",
        "model": model,
        "operation": operation,
        "points": points,
        "state": state,
        "rounds": 1,
        "time_samples_seconds": [median],
        "median_seconds": median,
        "fault_calls_total": calls,
        "fault_calls_top_level": 1,
        "field_calls_total": 3,
        "fault_points_total": fault_points,
        "field_points_total": 7,
        "peak_extra_rss_bytes": rss,
        "rss_unit": "bytes",
        "code_revision": code_revision,
        "dirty": False,
        "prerequisite_parent_revision": "parent",
        "prerequisite_fix_revision": "fix",
        "baseline_revision": "baseline",
        "input_hash": "0" * 64,
    }


def _document(rows: list[tuple], code_revision: str, **provenance_overrides) -> dict:
    """Return a complete schema v1 result document for the given rows."""
    return {
        "schema_version": 1,
        "cases": [_case(row, code_revision) for row in rows],
        "provenance": _provenance(**provenance_overrides),
        "pytest_benchmark": {"raw": {"stats": {"mean": 2.0}}},
    }


def _compliant_pair() -> tuple[dict, dict]:
    """Return a base and a new document that meet every target.

    The new document also has the extra inputs of the acceptance check:
    a repeat run (SC-007), the suite wall time (SC-008) and the gradient
    decision (OP1).
    """
    base = _document(BASE_ROWS, "baseline-rev")
    new = _document(NEW_ROWS, "new-rev")
    repeat = copy.deepcopy(new)
    for row in repeat["cases"]:
        row["median_seconds"] *= 1.05
        row["time_samples_seconds"] = [row["median_seconds"]]
    new["repeat"] = repeat
    new["suite_wall_seconds"] = 600.0
    new["gradient_decision"] = GRADIENT_RULE
    return base, new


def _row(document: dict, case_id: str) -> dict:
    """Return the row of one case id in a document."""
    rows = {row["case_id"]: row for row in document["cases"]}
    return rows[case_id]


def _set(document: dict, case_id: str, **fields) -> None:
    """Set fields of one case row in place."""
    _row(document, case_id).update(fields)


def test_compliant_pair_passes_every_target():
    """A pair that meets SC-002 to SC-010 passes, and each target has its limit."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    result = check_acceptance(base, new)
    assert result["passed"] is True
    assert set(result["criteria"]) == set(CRITERIA)
    for criterion in CRITERIA:
        row = result["criteria"][criterion]
        assert row["passed"] is True, criterion
        assert row["limit"] == pytest.approx(LIMITS[criterion]), criterion
        assert row["reason"] is None, criterion


@pytest.mark.parametrize(
    "criterion,mutate",
    [
        pytest.param(
            "SC-002",
            lambda base, new: _set(new, OLD_F4U3, median_seconds=1.01),
            id="SC-002-time-over-limit",
        ),
        pytest.param(
            "SC-002",
            lambda base, new: _set(new, OLD_F4U3, fault_calls_total=5),
            id="SC-002-calls-not-four",
        ),
        pytest.param(
            "SC-003",
            lambda base, new: _set(new, ALL_F4U3, median_seconds=1.3),
            id="SC-003-time-over-limit",
        ),
        pytest.param(
            "SC-003",
            lambda base, new: _set(new, ALL_F4U3, fault_calls_total=5),
            id="SC-003-calls-not-four",
        ),
        pytest.param(
            "SC-004",
            lambda base, new: _set(new, ALL_F4U0, median_seconds=0.6),
            id="SC-004-ratio-over-limit",
        ),
        pytest.param(
            "SC-005",
            lambda base, new: _set(new, PLUGIN_OLD, median_seconds=4.2),
            id="SC-005-time-over-limit",
        ),
        pytest.param(
            "SC-006",
            lambda base, new: _set(new, OLD_F4U3, peak_extra_rss_bytes=1_200_001),
            id="SC-006-memory-over-limit",
        ),
        pytest.param(
            "SC-007",
            lambda base, new: _set(new["repeat"], ALL_F4U3, median_seconds=1.2),
            id="SC-007-noise-over-limit",
        ),
        pytest.param(
            "SC-008",
            lambda base, new: new.update(suite_wall_seconds=900.0),
            id="SC-008-suite-not-under-900s",
        ),
        pytest.param(
            "SC-009",
            lambda base, new: _set(new, SINGLE_F1U0, median_seconds=0.35),
            id="SC-009-per-point-over-limit",
        ),
        pytest.param(
            "SC-010",
            lambda base, new: _set(new, WARM_F4U3, median_seconds=0.11),
            id="SC-010-warm-time-over-limit",
        ),
        pytest.param(
            "SC-010",
            lambda base, new: _set(new, WARM_F4U3, fault_calls_total=1),
            id="SC-010-warm-calls-not-zero",
        ),
    ],
)
def test_each_failed_target_fails_the_whole_result(criterion, mutate):
    """Each target that is not met fails its row and the overall result, with a reason."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    mutate(base, new)
    result = check_acceptance(base, new)
    assert result["criteria"][criterion]["passed"] is False
    assert result["criteria"][criterion]["reason"]
    assert result["passed"] is False


def test_target_limit_is_inclusive_for_sc_002():
    """A ratio exactly at the SC-002 limit (0.40) passes."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(new, OLD_F4U3, median_seconds=1.0)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-002"]["passed"] is True
    assert result["criteria"]["SC-002"]["value"] == pytest.approx(0.40)


def test_missing_case_is_an_explicit_failure():
    """A case missing from the new document fails the case check and its target."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new["cases"] = [row for row in new["cases"] if row["case_id"] != OLD_F4U3]
    result = check_acceptance(base, new)
    assert result["passed"] is False
    assert result["cases"]["passed"] is False
    assert any(OLD_F4U3 in text for text in result["cases"]["problems"])
    assert result["criteria"]["SC-002"]["passed"] is False


def test_missing_warm_case_fails_sc_010():
    """A missing warm case fails SC-010 explicitly, not by skipping it."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new["cases"] = [row for row in new["cases"] if row["case_id"] != WARM_F4U3]
    result = check_acceptance(base, new)
    assert result["passed"] is False
    assert result["criteria"]["SC-010"]["passed"] is False


def test_duplicate_case_in_new_fails_the_case_check():
    """A duplicate case id in the new document is a failure."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new["cases"].append(copy.deepcopy(_row(new, OLD_F4U3)))
    result = check_acceptance(base, new)
    assert result["passed"] is False
    assert result["cases"]["passed"] is False
    assert any("duplicate" in text.lower() for text in result["cases"]["problems"])


@pytest.mark.parametrize(
    "mutate",
    [
        pytest.param(
            lambda base, new: new["provenance"].update(python="3.12.4"),
            id="python-version-differs",
        ),
        pytest.param(
            lambda base, new: new.update(schema_version=2),
            id="schema-version-differs",
        ),
        pytest.param(
            lambda base, new: _set(new, OLD_F4U3, input_hash="1" * 64),
            id="input-hash-differs",
        ),
        pytest.param(
            lambda base, new: _set(new, OLD_F4U3, baseline_revision="other"),
            id="baseline-revision-differs",
        ),
    ],
)
def test_incompatible_provenance_fails(mutate):
    """Incompatible provenance or input identity fails the provenance check and the result."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    mutate(base, new)
    result = check_acceptance(base, new)
    assert result["provenance"]["passed"] is False
    assert result["provenance"]["problems"]
    assert result["passed"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        pytest.param(lambda base, new: new.pop("gradient_decision"), id="decision-missing"),
        pytest.param(
            lambda base, new: new.update(gradient_decision="one-restoration-per-call"),
            id="decision-other-rule",
        ),
    ],
)
def test_missing_gradient_decision_fails(mutate):
    """A missing or different gradient decision fails the result explicitly."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    mutate(base, new)
    result = check_acceptance(base, new)
    assert result["gradient_decision"]["passed"] is False
    assert result["gradient_decision"]["expected"] == GRADIENT_RULE
    assert result["passed"] is False


def test_missing_repeat_run_fails_sc_007():
    """Without a repeat run SC-007 cannot pass, and the result fails."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new.pop("repeat")
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is False
    assert result["passed"] is False


def test_missing_suite_time_fails_sc_008():
    """Without the suite wall time SC-008 cannot pass, and the result fails."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new.pop("suite_wall_seconds")
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-008"]["passed"] is False
    assert result["passed"] is False
