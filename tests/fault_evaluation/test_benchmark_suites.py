"""Tests for the standard and full benchmark suites and the pass case rules (D9, D10).

No test runs a real model or a real benchmark.
"""

from __future__ import annotations

import importlib

import pytest

run_module = importlib.import_module("benchmarks.fault_evaluation.run")
bench_module = importlib.import_module("benchmarks.fault_evaluation.bench_cases")

LARGE_MODELS = {
    "f0_u0",
    "f0_u1",
    "f0_u2",
    "f0_u3",
    "f1_u0",
    "f1_u1",
    "f1_u2",
    "f1_u3",
    "stopped_fault",
    "onlap_domain",
    "faulted_fold",
}
OPERATIONS = (
    "oldest_value",
    "all_group_values",
    "model_ids",
    "oldest_gradient",
    "model_gradient",
)


def _api(name: str):
    """Return an object of the runner, or fail with a clear reason when it is absent."""
    value = getattr(run_module, name, None)
    assert value is not None, f"benchmarks.fault_evaluation.run.{name} is required by D2"
    return value


def _ids(cases: list[dict]) -> list[str]:
    return [case["case_id"] for case in cases]


def test_full_suite_has_430_cold_cases():
    """enumerate_cases keeps its name and returns the full suite. Only the cold state remains."""
    cases = run_module.enumerate_cases()
    assert len(cases) == 430
    assert {case["state"] for case in cases} == {"cold"}
    assert {case["points"] for case in cases} == {5_000, 20_000, 50_000, 100_000, 500_000}


def test_standard_suite_counts():
    """The standard suite has 131 cases: 105 small, 10 extra and 16 large."""
    cases = _api("standard_cases")()
    assert len(cases) == 131
    by_size = {}
    for case in cases:
        by_size[case["points"]] = by_size.get(case["points"], 0) + 1
    assert by_size == {5_000: 105, 20_000: 10, 50_000: 9, 100_000: 7}
    assert len(_api("STANDARD_LARGE_CASES")) == 16


def test_standard_ids_are_in_full_list_in_same_order():
    """The standard list is a subsequence of the full list."""
    full_cases = run_module.enumerate_cases()
    full = _ids(full_cases)
    standard_cases = _api("standard_cases")()
    standard = _ids(standard_cases)
    assert len(set(standard)) == len(standard)
    positions = [full.index(case_id) for case_id in standard]
    assert positions == sorted(positions)
    by_id = {case["case_id"]: case for case in full_cases}
    assert all(by_id[case["case_id"]] == case for case in standard_cases)


def test_large_cases_use_only_the_allowed_models_and_literal_ids():
    """The 16 large ids use the 11 chosen models and are stored as literal text."""
    large = list(_api("STANDARD_LARGE_CASES"))
    assert all(isinstance(case_id, str) for case_id in large)
    full = {case["case_id"]: case for case in run_module.enumerate_cases()}
    assert set(large) <= set(full)
    assert {full[case_id]["model"] for case_id in large} <= LARGE_MODELS
    assert {full[case_id]["points"] for case_id in large} == {50_000, 100_000}
    standard_large = [
        case["case_id"] for case in _api("standard_cases")() if case["points"] in (50_000, 100_000)
    ]
    assert sorted(standard_large) == sorted(large)


def test_standard_suite_contains_the_acceptance_cases():
    """Every case that SC-002 to SC-010 use is in the standard suite."""
    standard = set(_ids(_api("standard_cases")()))
    targets = (("f4_u3", 5_000), ("f4_u0", 5_000), ("plugin", 20_000), ("f1_u0", 20_000))
    wanted = {
        f"{model}:{operation}:{points}:{state}"
        for model, points in targets
        for operation in OPERATIONS
        for state in ("cold",)
    }
    assert wanted <= standard


def test_select_cases_names():
    """select_cases picks the standard or full list and rejects other names."""
    select_cases = _api("select_cases")
    assert _ids(select_cases("standard")) == _ids(_api("standard_cases")())
    assert _ids(select_cases("full")) == _ids(run_module.enumerate_cases())
    with pytest.raises(ValueError):
        select_cases("other")


EXPECTED_LARGE_CASES = (
    "f0_u0:all_group_values:100000:cold",
    "f0_u1:oldest_value:100000:cold",
    "f0_u1:model_ids:50000:cold",
    "f0_u2:oldest_value:50000:cold",
    "f0_u3:model_ids:100000:cold",
    "f1_u0:all_group_values:100000:cold",
    "f1_u1:oldest_value:50000:cold",
    "f1_u2:oldest_value:100000:cold",
    "f1_u2:oldest_gradient:50000:cold",
    "f1_u3:oldest_value:50000:cold",
    "f1_u3:all_group_values:100000:cold",
    "stopped_fault:model_ids:50000:cold",
    "stopped_fault:oldest_gradient:50000:cold",
    "onlap_domain:model_ids:50000:cold",
    "faulted_fold:oldest_value:50000:cold",
    "faulted_fold:all_group_values:100000:cold",
)


def test_states_are_only_cold():
    """The benchmark has no warm state (D9)."""
    assert _api("STATES") == ("cold",)


def test_large_cases_are_the_16_ids_of_plan_d9_in_order():
    """STANDARD_LARGE_CASES is the literal list of plan D9."""
    assert tuple(_api("STANDARD_LARGE_CASES")) == EXPECTED_LARGE_CASES


def test_standard_suite_has_no_onlap_gradient_case_at_100000_points():
    """The five gradient cases at 100,000 points are not in the standard list."""
    standard = _ids(_api("standard_cases")())
    assert not [
        case_id
        for case_id in standard
        if case_id.startswith("onlap_domain:") and ":100000:" in case_id and "gradient" in case_id
    ]


def test_pass_points_table():
    """Timing runs all cases. Counters run 5,000 and 20,000 points. Memory runs 5,000."""
    assert _api("PASS_POINTS") == {"timing": None, "counters": (5_000, 20_000), "rss": (5_000,)}


def _sized_cases() -> list[dict]:
    return [
        {"case_id": f"m:oldest_value:{points}:cold", "points": points}
        for points in (5_000, 20_000, 50_000, 100_000, 500_000)
    ]


def test_pass_cases_full_gives_all_cases_for_every_pass():
    """The full suite runs each pass on every case."""
    pass_cases = _api("pass_cases")
    cases = _sized_cases()
    for pass_kind in ("timing", "counters", "rss"):
        assert pass_cases(cases, pass_kind, "full") == cases


def test_pass_cases_standard_filters_by_points():
    """The standard suite filters the counter and memory passes by point count."""
    pass_cases = _api("pass_cases")
    cases = _sized_cases()
    assert pass_cases(cases, "timing", "standard") == cases
    assert [case["points"] for case in pass_cases(cases, "counters", "standard")] == [5_000, 20_000]
    assert [case["points"] for case in pass_cases(cases, "rss", "standard")] == [5_000]


def test_pass_cases_rejects_unknown_pass_and_suite():
    """An unknown pass or suite name raises ValueError."""
    pass_cases = _api("pass_cases")
    with pytest.raises(ValueError):
        pass_cases(_sized_cases(), "other", "standard")
    with pytest.raises(ValueError):
        pass_cases(_sized_cases(), "timing", "other")


def _first_pass_ids(monkeypatch, tmp_path, extra_argv):
    """Run main with a fake pass that stops the run. Return the case ids of the first pass."""
    full = [
        {
            "case_id": "tiny:full",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2,
            "state": "cold",
        }
    ]
    # main must call enumerate_cases through the module global for the full suite,
    # so that this fake list is what "--suite full" runs.
    monkeypatch.setattr(run_module, "enumerate_cases", lambda: full)
    seen = []

    def fake_pass(pass_kind, selected_cases, *, output_path, rounds):
        seen.append(_ids(selected_cases))
        raise RuntimeError("stop after the first pass")

    monkeypatch.setattr(run_module, "_run_fresh_process_pass", fake_pass)
    argv = ["--output", str(tmp_path / "out.json"), "--work-dir", str(tmp_path / "work")]
    with pytest.raises(RuntimeError, match="stop after the first pass"):
        run_module.main([*argv, *extra_argv])
    return seen[0]


def test_main_suite_option_defaults_to_standard(monkeypatch, tmp_path):
    """Without --suite, main runs the standard list."""
    ids = _first_pass_ids(monkeypatch, tmp_path, [])
    assert ids == _ids(_api("standard_cases")())


def test_main_suite_full_runs_the_full_list(monkeypatch, tmp_path):
    """--suite full runs what enumerate_cases returns."""
    assert _first_pass_ids(monkeypatch, tmp_path, ["--suite", "full"]) == ["tiny:full"]


def test_main_rejects_unknown_suite(tmp_path):
    """An unknown suite name is a usage error."""
    with pytest.raises(SystemExit):
        run_module.main(["--output", str(tmp_path / "out.json"), "--suite", "other"])


class _FakeBenchmark:
    """Replace pytest-benchmark: run setup, then the call, once for each round."""

    def __init__(self):
        self.extra_info = {}

    def pedantic(self, call, *, setup, rounds, warmup_rounds):
        assert warmup_rounds == 0
        for _ in range(rounds):
            setup()
            call()


class _NoCacheRead:
    """A model that fails the test if run_case reads its evaluation cache."""

    @property
    def evaluation_cache(self):
        raise AssertionError("run_case must not read evaluation_cache (D9)")


def test_run_case_never_primes_and_never_reads_the_cache(monkeypatch):
    """A case makes one timed call in each round. It never clears or primes a cache."""
    events = []
    monkeypatch.setitem(
        bench_module._OPERATIONS, "oldest_value", lambda model, points: events.append("call")
    )
    for state in ("cold", "warm"):
        events.clear()
        case = {"operation": "oldest_value", "state": state}
        bench_module.run_case(_FakeBenchmark(), _NoCacheRead(), None, case, rounds=3)
        assert events == ["call"] * 3, state
