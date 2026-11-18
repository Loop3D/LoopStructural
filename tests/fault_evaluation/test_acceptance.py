"""Contract tests for numerical performance acceptance reporting (schema v1).

The acceptance check reads a base and a new schema v1 result document. It
reports one row for each target SC-002 to SC-009. A failed target is a
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
}

# Case ids of the cases that the targets use. The id format matches run.py.
OLD_F4U3 = "f4_u3:oldest_value:5000:cold"
ALL_F4U3 = "f4_u3:all_group_values:5000:cold"
ALL_F4U0 = "f4_u0:all_group_values:5000:cold"
PLUGIN_OLD = "plugin:oldest_value:20000:cold"
SINGLE_F1U0 = "f1_u0:oldest_value:20000:cold"

# Row fields: model, operation, points, state, median seconds, fault calls,
# fault points, peak extra RSS bytes.
BASE_ROWS = [
    ("f4_u3", "oldest_value", 5000, "cold", 2.5, 16, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "cold", 8.65, 52, 5000, 1_000_000),
    ("f4_u0", "all_group_values", 5000, "cold", 0.67, 12, 5000, 1_000_000),
    ("plugin", "oldest_value", 20000, "cold", 8.2, 6, 20000, 1_000_000),
    ("f1_u0", "oldest_value", 20000, "cold", 0.68, 2, 20000, 500_000),
]
NEW_ROWS = [
    ("f4_u3", "oldest_value", 5000, "cold", 1.0, 4, 5000, 1_000_000),
    ("f4_u3", "all_group_values", 5000, "cold", 1.0, 4, 5000, 1_000_000),
    ("f4_u0", "all_group_values", 5000, "cold", 0.8, 4, 5000, 1_000_000),
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
    """A pair that meets SC-002 to SC-009 passes, and each target has its limit."""
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
            lambda base, new: (
                _set(base, OLD_F4U3, peak_extra_rss_bytes=100_000_000),
                _set(new, OLD_F4U3, peak_extra_rss_bytes=120_000_001),
            ),
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


class _Fault:
    def __init__(self, name):
        self.name = name


class _Feature:
    def __init__(self, *fault_names):
        self.faults = [_Fault(name) for name in fault_names]


class _Model:
    def __init__(self, features):
        self.features = features


def test_distinct_fault_sets_counts_distinct_ordered_tuples():
    """Features with equal ordered fault names count once. A reversed order counts twice."""
    count = _api(ACCEPTANCE_MODULE, "_distinct_fault_sets")
    model = _Model(
        [
            _Feature("a", "b"),
            _Feature("a", "b"),
            _Feature("b", "a"),
            _Feature(),
            _Feature("c"),
        ]
    )
    assert count(model) == 3


def test_distinct_fault_sets_is_zero_without_faults():
    """A model with no faults has 0 fault sets."""
    count = _api(ACCEPTANCE_MODULE, "_distinct_fault_sets")
    assert count(_Model([_Feature(), _Feature()])) == 0
    assert count(_Model([])) == 0


def test_memory_accounting_fields():
    """The FR-010 accounting uses the distinct fault sets and the restored batches."""
    accounting = _api(ACCEPTANCE_MODULE, "_fr010_accounting")
    fields = accounting(points=100, fault_sets=2, batches={100: 2, 50: 1}, peak=9600)
    assert fields["distinct_fault_sets"] == 2
    assert fields["allowed_bytes"] == 24 * 100 * 2
    assert fields["restored_bytes_observed"] == 24 * (200 + 50)
    assert fields["peak_over_allowed"] == pytest.approx(9600 / 4800)
    assert accounting(points=100, fault_sets=0, batches={}, peak=5)["peak_over_allowed"] is None


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


def _short_rows(result: dict) -> dict:
    """Return the short_cases rows of a result, keyed by case id."""
    assert "short_cases" in result, "result['short_cases'] is required by D2b"
    return {row["case_id"]: row for row in result["short_cases"]}


def test_sc_007_checks_a_long_case_and_names_it():
    """A case with a baseline median of 1.0 s or more and noise 1.10 or more fails SC-007."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(new["repeat"], PLUGIN_OLD, median_seconds=3.0 * 1.10)
    result = check_acceptance(base, new)
    row = result["criteria"]["SC-007"]
    assert row["passed"] is False
    assert PLUGIN_OLD in row["reason"]
    assert result["passed"] is False


def test_sc_007_limit_is_a_strict_bound_for_a_long_case():
    """A noise ratio just below 1.10 passes SC-007 for a case of 1.0 s or more."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(new["repeat"], PLUGIN_OLD, median_seconds=3.0 * 1.09)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is True
    assert PLUGIN_OLD not in _short_rows(result)


def test_sc_007_short_case_is_listed_and_does_not_fail():
    """A case with a baseline median below 1.0 s and noise 1.5 is listed, not failed."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    assert _row(base, ALL_F4U0)["median_seconds"] < 1.0
    _set(new["repeat"], ALL_F4U0, median_seconds=0.8 * 1.5)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is True
    assert result["criteria"]["SC-007"]["reason"] is None
    rows = _short_rows(result)
    assert ALL_F4U0 in rows
    assert rows[ALL_F4U0]["noise_ratio"] == pytest.approx(1.5)


def test_sc_007_short_cases_list_every_short_case_and_no_long_case():
    """The short_cases list holds each case with a baseline median below 1.0 s."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    result = check_acceptance(base, new)
    rows = _short_rows(result)
    assert set(rows) == {ALL_F4U0, SINGLE_F1U0}
    for row in rows.values():
        assert row["noise_ratio"] == pytest.approx(1.05)


def test_sc_007_baseline_median_of_1_0_s_is_checked_not_short():
    """A baseline median of exactly 1.0 s is a checked case, not a short case."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(base, ALL_F4U0, median_seconds=1.0)
    _set(new["repeat"], ALL_F4U0, median_seconds=0.8 * 1.5)
    result = check_acceptance(base, new)
    assert ALL_F4U0 not in _short_rows(result)
    assert result["criteria"]["SC-007"]["passed"] is False
    assert ALL_F4U0 in result["criteria"]["SC-007"]["reason"]


def test_sc_007_missing_repeat_case_still_fails():
    """A case that is missing from the repeat document fails SC-007."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new["repeat"]["cases"] = [row for row in new["repeat"]["cases"] if row["case_id"] != PLUGIN_OLD]
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is False
    assert PLUGIN_OLD in result["criteria"]["SC-007"]["reason"]
    assert result["passed"] is False


def test_sc_007_zero_median_in_a_long_case_fails():
    """A zero median in new or repeat for a case of 1.0 s or more fails SC-007.

    A zero median is an invalid measurement for every case, long or short.
    """
    check_acceptance = _check_acceptance()
    for document in ("new", "repeat"):
        base, new = _compliant_pair()
        target = new if document == "new" else new["repeat"]
        _set(target, PLUGIN_OLD, median_seconds=0.0)
        result = check_acceptance(base, new)
        assert result["criteria"]["SC-007"]["passed"] is False, document
        assert result["passed"] is False, document


def test_sc_007_zero_baseline_median_fails_and_names_the_case():
    """A zero baseline median is an invalid measurement. It fails SC-007 and names the case.

    This holds for a case that would otherwise be short.
    """
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(base, SINGLE_F1U0, median_seconds=0.0)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is False
    assert SINGLE_F1U0 in result["criteria"]["SC-007"]["reason"]
    assert result["passed"] is False


def test_sc_007_zero_repeat_median_in_a_short_case_fails_and_names_the_case():
    """A zero median in the repeat document fails SC-007 for a short case and names it."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    assert _row(base, ALL_F4U0)["median_seconds"] < 1.0
    _set(new["repeat"], ALL_F4U0, median_seconds=0.0)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-007"]["passed"] is False
    assert ALL_F4U0 in result["criteria"]["SC-007"]["reason"]
    assert result["passed"] is False


def test_missing_suite_time_fails_sc_008():
    """Without the suite wall time SC-008 cannot pass, and the result fails."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    new.pop("suite_wall_seconds")
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-008"]["passed"] is False
    assert result["passed"] is False


def test_there_is_no_sc_010_row():
    """SC-010 is removed with the warm cases (D10)."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    assert "SC-010" not in check_acceptance(base, new)["criteria"]


def test_memory_modes_and_allowance():
    """The large-memory check has only the cold mode. The allowance is 16 MiB (D10)."""
    assert _api(ACCEPTANCE_MODULE, "MEMORY_MODES") == ("cold",)
    assert _api(ACCEPTANCE_MODULE, "MEMORY_ALLOWANCE_BYTES") == 16 * 1024 * 1024


def _sc006(base_bytes, new_bytes) -> dict:
    """Return the SC-006 row for one baseline and one candidate memory value."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(base, OLD_F4U3, peak_extra_rss_bytes=base_bytes)
    _set(new, OLD_F4U3, peak_extra_rss_bytes=new_bytes)
    return check_acceptance(base, new)["criteria"]["SC-006"]


def test_sc_006_passes_a_small_baseline_within_the_allowance():
    """1.5 x a 1 MB baseline is over 120 %, but the extra 0.5 MB is under the allowance."""
    assert _sc006(1_000_000, 1_500_000)["passed"] is True


def test_sc_006_fails_a_large_baseline_over_both_limits():
    """1.5 x a 100 MB baseline is over 120 % and over the allowance."""
    row = _sc006(100_000_000, 150_000_000)
    assert row["passed"] is False
    assert row["reason"]


def test_sc_006_passes_at_the_ratio_limit_for_a_large_baseline():
    """A ratio of 1.20 passes for a large baseline."""
    assert _sc006(100_000_000, 120_000_000)["passed"] is True


def test_sc_006_allowance_limit_is_inclusive_and_a_zero_baseline_uses_only_it():
    """The extra memory of exactly the allowance passes. One byte more fails."""
    allowance = _api(ACCEPTANCE_MODULE, "MEMORY_ALLOWANCE_BYTES")
    assert _sc006(0, allowance)["passed"] is True
    assert _sc006(0, allowance + 1)["passed"] is False
    assert _sc006(100_000_000, 100_000_000 + allowance)["passed"] is True
    assert _sc006(100_000_000, 100_000_000 + allowance + 1)["passed"] is False


def test_candidate_without_memory_value_is_not_a_failure():
    """A candidate case that the memory pass did not run (None) does not fail SC-006."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(new, SINGLE_F1U0, peak_extra_rss_bytes=None)
    result = check_acceptance(base, new)
    assert result["criteria"]["SC-006"]["passed"] is True
    assert result["passed"] is True


@pytest.mark.parametrize(
    "case_id,field,pass_name",
    [
        pytest.param(OLD_F4U3, "fault_calls_total", "counters", id="SC-002-calls"),
        pytest.param(ALL_F4U3, "fault_calls_total", "counters", id="SC-003-calls"),
        pytest.param(SINGLE_F1U0, "fault_points_total", "counters", id="SC-009-points"),
    ],
)
def test_criterion_that_needs_a_missing_counter_fails_with_the_pass_name(case_id, field, pass_name):
    """A needed counter of None fails the criterion with the case and the pass."""
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    _set(new, case_id, **{field: None})
    result = check_acceptance(base, new)
    assert result["passed"] is False
    reasons = [row["reason"] for row in result["criteria"].values() if row["reason"]]
    assert any(
        case_id in reason and f"was not measured in the {pass_name} pass" in reason
        for reason in reasons
    ), reasons


def test_baseline_with_warm_ids_goes_through_cold_view():
    """A baseline of the code before this feature (C-016) has warm rows.

    A warm row with a cold row of the same case goes. A lone warm row becomes the
    cold row.
    """
    check_acceptance = _check_acceptance()
    base, new = _compliant_pair()
    plugin = _row(base, PLUGIN_OLD)
    plugin.update(case_id=PLUGIN_OLD.replace(":cold", ":warm"), state="warm")
    extra_warm = copy.deepcopy(_row(base, ALL_F4U3))
    extra_warm.update(case_id=ALL_F4U3.replace(":cold", ":warm"), state="warm", median_seconds=99.0)
    base["cases"].append(extra_warm)
    result = check_acceptance(base, new)
    assert result["cases"]["problems"] == []
    assert result["passed"] is True
