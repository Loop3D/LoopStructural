"""Tests for the standard and full benchmark suites and the warm priming rule.

No test runs a real model or a real benchmark.
"""

from __future__ import annotations

import importlib
from types import SimpleNamespace

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


def test_full_suite_has_860_cases():
    """enumerate_cases keeps its name and returns the full suite."""
    cases = run_module.enumerate_cases()
    assert len(cases) == 860
    assert {case["points"] for case in cases} == {5_000, 20_000, 50_000, 100_000, 500_000}


def test_standard_suite_counts():
    """The standard suite has 254 cases: 210 small, 20 extra and 24 large."""
    cases = _api("standard_cases")()
    assert len(cases) == 254
    by_size = {}
    for case in cases:
        by_size[case["points"]] = by_size.get(case["points"], 0) + 1
    assert by_size == {5_000: 210, 20_000: 20, 50_000: 12, 100_000: 12}
    assert len(_api("STANDARD_LARGE_CASES")) == 24


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
    """The 24 large ids use the 11 chosen models and are stored as literal text."""
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
        for state in ("cold", "warm")
    }
    assert wanted <= standard


def test_select_cases_names():
    """select_cases picks the standard or full list and rejects other names."""
    select_cases = _api("select_cases")
    assert _ids(select_cases("standard")) == _ids(_api("standard_cases")())
    assert _ids(select_cases("full")) == _ids(run_module.enumerate_cases())
    with pytest.raises(ValueError):
        select_cases("other")


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


class _FakeCache:
    def __init__(self, events):
        self.events = events

    def clear(self):
        self.events.append("clear")


def _run_fake_case(monkeypatch, state, *, with_cache, rounds=3):
    """Run one fake case. Return the list of cache clears and operation calls."""
    events = []
    model = SimpleNamespace()
    if with_cache:
        model.evaluation_cache = _FakeCache(events)
    monkeypatch.setitem(
        bench_module._OPERATIONS, "oldest_value", lambda model, points: events.append("call")
    )
    case = {"operation": "oldest_value", "state": state}
    bench_module.run_case(_FakeBenchmark(), model, None, case, rounds=rounds)
    return events


def test_warm_case_without_cache_makes_no_priming_call(monkeypatch):
    """A model with no evaluation_cache has nothing to fill: one call in each round."""
    assert _run_fake_case(monkeypatch, "warm", with_cache=False) == ["call"] * 3


def test_warm_case_with_cache_primes_after_each_clear(monkeypatch):
    """A model with a cache is cleared and primed once in each round, then timed."""
    events = _run_fake_case(monkeypatch, "warm", with_cache=True)
    assert events == ["clear", "call", "call"] * 3


def test_cold_case_clears_the_cache_and_never_primes(monkeypatch):
    """A cold case clears the cache in each round and makes only the timed call."""
    assert _run_fake_case(monkeypatch, "cold", with_cache=True) == ["clear", "call"] * 3
    assert _run_fake_case(monkeypatch, "cold", with_cache=False) == ["call"] * 3
