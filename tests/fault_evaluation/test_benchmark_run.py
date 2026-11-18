"""Contract tests for opt-in fault evaluation benchmark tooling."""

from __future__ import annotations

import copy
import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


def _api(module_name: str, name: str):
    """Import a focused API and report an actionable absent-tool failure."""
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as error:
        if error.name == module_name or module_name.startswith(f"{error.name}."):
            pytest.fail(f"{module_name}.{name} is required by D2; benchmark tools are missing")
        raise
    value = getattr(module, name, None)
    assert value is not None, f"{module_name}.{name} is required by D2"
    return value


def _tiny_result(pass_kind: str, cases: list[dict]) -> dict:
    """Build a valid two-case result for a mocked isolated runner pass."""
    timing = pass_kind == "timing"
    counter = pass_kind == "counters"
    rows = []
    for case in cases:
        row = {
            **case,
            "rounds": 2,
            "time_samples_seconds": [0.01, 0.02] if timing else [],
            "median_seconds": 0.015 if timing else None,
            "fault_calls_total": 2 if counter else 0,
            "fault_calls_top_level": 1 if counter else 0,
            "field_calls_total": 3 if counter else 0,
            "fault_points_total": 4 if counter else 0,
            "field_points_total": 7 if counter else 0,
            "peak_extra_rss_bytes": 4096 if pass_kind == "rss" else 0,
            "rss_unit": "bytes",
            "code_revision": "abc123",
            "dirty": False,
            "prerequisite_parent_revision": "parent",
            "prerequisite_fix_revision": "fix",
            "baseline_revision": "baseline",
            "input_hash": "0" * 64,
        }
        rows.append(row)
    return {
        "schema_version": 1,
        "pass": pass_kind,
        "cases": rows,
        "provenance": {
            "python": "3.x",
            "numpy": "2.x",
            "scipy": "1.x",
            "loopstructural": "2.x",
            "os": "test",
            "cpu": "test cpu",
            "ram_bytes": 1_000_000,
        },
        "pytest_benchmark": {"raw": {"stats": {"mean": 0.015}}} if timing else {},
    }


def test_enumeration_is_stable_and_covers_full_matrix_and_addons():
    """All required models/operations/sizes have stable cold cases (D9)."""
    enumerate_cases = _api("benchmarks.fault_evaluation.run", "enumerate_cases")
    cases = enumerate_cases()
    repeated = enumerate_cases()
    assert [item["case_id"] for item in cases] == [item["case_id"] for item in repeated]
    assert len({item["case_id"] for item in cases}) == len(cases)

    models = {f"f{faults}_u{units}" for faults in (0, 1, 2, 4) for units in range(4)}
    models |= {"stopped_fault", "cutting_faults", "onlap_domain", "faulted_fold", "plugin"}
    operations = {
        "oldest_value",
        "all_group_values",
        "model_ids",
        "oldest_gradient",
        "model_gradient",
    }
    sizes = {5_000, 50_000, 500_000}
    required = {
        (model, operation, count, state)
        for model in models
        for operation in operations
        for count in sizes
        for state in ("cold",)
    }
    observed = {(item["model"], item["operation"], item["points"], item["state"]) for item in cases}
    assert required <= observed
    assert {item["model"] for item in cases} == models

    for operation in operations:
        assert ("plugin", operation, 20_000, "cold") in observed
        assert ("f1_u0", operation, 20_000, "cold") in observed
    assert {item["state"] for item in cases} == {"cold"}
    assert len(cases) == 430


def test_instrumentation_counts_nested_faults_and_local_fields_and_restores(monkeypatch):
    """Only opt-in fault/local-field wrappers count; every patch restores on error."""
    from LoopStructural.modelling.features.fault import FaultSegment
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator

    # Replace implementations with tiny bound-call-safe fakes before installing
    # instrumentation. The production evaluation methods are never run.
    monkeypatch.setattr(
        FaultSegment,
        "apply_to_points",
        lambda self, points, *args, **kwargs: points,
    )
    monkeypatch.setattr(
        DiscreteInterpolator,
        "_evaluate_value_local",
        lambda self, points: np.zeros(len(points)),
    )
    monkeypatch.setattr(
        DiscreteInterpolator,
        "_evaluate_gradient_local",
        lambda self, points: np.zeros((len(points), 3)),
    )
    instrumentation_type = _api(
        "benchmarks.fault_evaluation.instrumentation", "EvaluationInstrumentation"
    )
    original_fault = FaultSegment.apply_to_points
    original_value = DiscreteInterpolator._evaluate_value_local
    original_gradient = DiscreteInterpolator._evaluate_gradient_local

    class FakeInterpolator:
        pass

    class FakeFault:
        pass

    class ChildFault(FakeFault):
        pass

    fault = FakeFault()
    child = ChildFault()
    points = np.zeros((3, 3))

    def nested_apply(self, values, *args, **kwargs):
        if self is fault:
            FaultSegment.apply_to_points(child, values[:1])
        return values

    monkeypatch.setattr(FaultSegment, "apply_to_points", nested_apply)
    instance = instrumentation_type()
    with pytest.raises(RuntimeError, match="instrumentation sentinel"):
        with instance:
            FaultSegment.apply_to_points(fault, points)
            DiscreteInterpolator._evaluate_value_local(FakeInterpolator(), points)
            DiscreteInterpolator._evaluate_gradient_local(FakeInterpolator(), np.zeros((4, 3)))
            raise RuntimeError("instrumentation sentinel")

    assert FaultSegment.apply_to_points is nested_apply
    assert DiscreteInterpolator._evaluate_value_local is original_value
    assert DiscreteInterpolator._evaluate_gradient_local is original_gradient
    counts = instance.snapshot()
    assert counts["fault_calls_total"] == 2
    assert counts["fault_calls_top_level"] == 1
    assert counts["fault_points_total"] == 4
    assert counts["field_calls_total"] == 2
    assert counts["field_points_total"] == 7


def test_result_validation_rejects_missing_provenance_counters_and_bad_ids():
    """Validation rejects incomplete schema, duplicate cases and invalid RSS units."""
    validate_result = _api("benchmarks.fault_evaluation.run", "validate_result")
    case = {
        "case_id": "tiny:oldest_value:2:cold",
        "model": "tiny",
        "operation": "oldest_value",
        "points": 2,
        "state": "cold",
        "rounds": 2,
        "time_samples_seconds": [0.01, 0.02],
        "median_seconds": 0.015,
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
    result = {
        "schema_version": 1,
        "cases": [case],
        "provenance": {
            "python": "3.x",
            "numpy": "2.x",
            "scipy": "1.x",
            "loopstructural": "2.x",
            "os": "test",
            "cpu": "test cpu",
            "ram_bytes": 1_000_000,
        },
        "pytest_benchmark": {"raw": {"stats": {"mean": 0.015}}},
    }
    validate_result(result)

    missing = json.loads(json.dumps(result))
    del missing["cases"][0]["baseline_revision"]
    with pytest.raises((ValueError, KeyError), match="baseline_revision|revision|field"):
        validate_result(missing)
    missing = json.loads(json.dumps(result))
    del missing["provenance"]["cpu"]
    with pytest.raises((ValueError, KeyError), match="cpu|provenance|field"):
        validate_result(missing)
    missing = json.loads(json.dumps(result))
    del missing["cases"][0]["field_points_total"]
    with pytest.raises((ValueError, KeyError), match="field_points_total|counter|field"):
        validate_result(missing)
    duplicate = json.loads(json.dumps(result))
    duplicate["cases"].append(dict(case))
    with pytest.raises(ValueError, match="duplicate|unique|case"):
        validate_result(duplicate)
    bad_units = json.loads(json.dumps(result))
    bad_units["cases"][0]["rss_unit"] = "pages"
    with pytest.raises(ValueError, match="rss|unit|bytes"):
        validate_result(bad_units)


def test_normal_benchmark_import_does_not_patch_runtime_methods():
    """Importing the runner alone leaves fault and interpolator methods intact."""
    from LoopStructural.modelling.features.fault import FaultSegment
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator

    original_fault = FaultSegment.apply_to_points
    original_value = DiscreteInterpolator._evaluate_value_local
    original_gradient = DiscreteInterpolator._evaluate_gradient_local
    _api("benchmarks.fault_evaluation.run", "main")
    assert FaultSegment.apply_to_points is original_fault
    assert DiscreteInterpolator._evaluate_value_local is original_value
    assert DiscreteInterpolator._evaluate_gradient_local is original_gradient


def test_timing_pass_uses_only_benchmark_suite_and_parses_raw_json(monkeypatch, tmp_path):
    """Timing delegates to the explicit suite launcher and parses its raw output."""
    run_pass = _api("benchmarks.fault_evaluation.run", "_run_fresh_process_pass")
    run_module = importlib.import_module(run_pass.__module__)
    cases = [
        {
            "case_id": "tiny:cold",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2,
            "state": "cold",
        }
    ]
    output = tmp_path / "timing.json"
    launches = []

    def fake_launch(selected_suite, args):
        launches.append((selected_suite, args))
        benchmark_path = Path(args[args.index("--benchmark-json") + 1])
        benchmark_path.write_text(
            json.dumps(
                {
                    "benchmarks": [
                        {
                            "extra_info": {"case_id": "tiny:cold"},
                            "stats": {
                                "data": [0.01, 0.02],
                                "median": 0.015,
                                "rounds": 2,
                                "iterations": 1,
                            },
                        }
                    ]
                }
            ),
            encoding="utf8",
        )
        return 0

    monkeypatch.setattr(run_module, "_launch_pytest", fake_launch)
    payload = run_pass("timing", cases, output_path=output, rounds=2)
    assert launches and Path(launches[0][0]).name == "bench_cases.py"
    assert payload["cases"][0]["time_samples_seconds"] == [0.01, 0.02]


def test_launcher_builds_fresh_python_pytest_command_without_child_process(monkeypatch):
    """Launcher command explicitly invokes this interpreter and suite; no child starts."""
    run_launcher = _api("benchmarks.fault_evaluation.run", "_launch_pytest")
    run_module = importlib.import_module(run_launcher.__module__)
    suite = Path("benchmarks/fault_evaluation/bench_cases.py").resolve()
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0)

    def forbidden_popen(*args, **kwargs):
        pytest.fail("launcher attempted an unmocked child process")

    monkeypatch.setattr(run_module.subprocess, "run", fake_run)
    monkeypatch.setattr(run_module.subprocess, "Popen", forbidden_popen)
    assert run_launcher(suite, ["--benchmark-json", "raw.json"]) == 0
    assert commands
    command = commands[0]
    assert command[:3] == [sys.executable, "-m", "pytest"]
    assert str(suite) in command
    assert "--benchmark-json" in command and "raw.json" in command


def _tiny_cases(*points: int) -> list[dict]:
    """Return one cold case for each point count."""
    return [
        {
            "case_id": f"tiny:oldest_value:{count}:cold",
            "model": "tiny",
            "operation": "oldest_value",
            "points": count,
            "state": "cold",
        }
        for count in points
    ]


def _run_main_with_fake_passes(monkeypatch, tmp_path, argv, cases, *, patch="enumerate_cases"):
    """Run main with fake passes. Return the calls and the merged document."""
    run_main = _api("benchmarks.fault_evaluation.run", "main")
    run_module = importlib.import_module(run_main.__module__)
    monkeypatch.setattr(run_module, patch, lambda: cases)
    calls = []

    def fake_pass(pass_kind, selected_cases, *, output_path, rounds):
        calls.append((pass_kind, [item["case_id"] for item in selected_cases], rounds))
        payload = _tiny_result(pass_kind, selected_cases)
        Path(output_path).write_text(json.dumps(payload), encoding="utf8")
        return payload

    monkeypatch.setattr(run_module, "_run_fresh_process_pass", fake_pass)
    monkeypatch.setattr(run_module, "_launch_pytest", lambda *args: 0)
    output = tmp_path / "merged.json"
    assert run_main(["--output", str(output), *argv]) == 0
    return calls, json.loads(output.read_text(encoding="utf8"))


def test_runner_separates_passes_and_writes_merged_case_metrics(monkeypatch, tmp_path):
    """Three mocked passes merge into schema output with timing raw statistics."""
    cases = _tiny_cases(2, 3)
    calls, merged = _run_main_with_fake_passes(
        monkeypatch, tmp_path, ["--rounds", "2", "--suite", "full"], cases
    )
    ids = [case["case_id"] for case in cases]
    assert [call[0] for call in calls] == ["timing", "counters", "rss"]
    assert all(call[1] == ids and call[2] == 2 for call in calls)
    assert merged["schema_version"] == 1
    assert merged["suite"] == "full"
    assert [row["case_id"] for row in merged["cases"]] == ids
    assert all(row["fault_calls_total"] == 2 for row in merged["cases"])
    assert all(row["field_points_total"] == 7 for row in merged["cases"])
    assert all(row["peak_extra_rss_bytes"] == 4096 for row in merged["cases"])
    assert merged["pytest_benchmark"]["raw"]["stats"]["mean"] == 0.015


def test_standard_run_gives_each_pass_its_own_cases_and_merges_missing_as_none(
    monkeypatch, tmp_path
):
    """The counter and memory passes run fewer cases. The merge writes None for the rest (D10)."""
    cases = _tiny_cases(5_000, 20_000, 50_000, 100_000)
    calls, merged = _run_main_with_fake_passes(
        monkeypatch, tmp_path, ["--suite", "standard"], cases, patch="standard_cases"
    )
    ids = [case["case_id"] for case in cases]
    assert [call[:2] for call in calls] == [
        ("timing", ids),
        ("counters", ids[:2]),
        ("rss", ids[:1]),
    ]
    assert merged["suite"] == "standard"
    rows = {row["case_id"]: row for row in merged["cases"]}
    assert list(rows) == ids
    counter_fields = _api("benchmarks.fault_evaluation.run", "COUNTER_FIELDS")
    for case_id in ids[2:]:
        assert all(rows[case_id][field] is None for field in counter_fields)
    for case_id in ids[:2]:
        assert all(rows[case_id][field] is not None for field in counter_fields)
    for case_id in ids[1:]:
        assert rows[case_id]["peak_extra_rss_bytes"] is None
        assert rows[case_id]["rss_unit"] == "bytes"
    assert rows[ids[0]]["peak_extra_rss_bytes"] == 4096


def _valid_cold_result() -> dict:
    """Return a valid result with one cold case."""
    tiny = _tiny_result("timing", _tiny_cases(2))
    tiny.pop("pass")
    return tiny


def test_result_validation_accepts_none_only_for_counters_and_memory():
    """A case that a pass did not measure has None counters or None memory (D10)."""
    validate_result = _api("benchmarks.fault_evaluation.run", "validate_result")
    counter_fields = _api("benchmarks.fault_evaluation.run", "COUNTER_FIELDS")
    result = _valid_cold_result()
    validate_result(result)
    for field in counter_fields:
        changed = copy.deepcopy(result)
        changed["cases"][0][field] = None
        validate_result(changed)
    changed = copy.deepcopy(result)
    changed["cases"][0]["peak_extra_rss_bytes"] = None
    validate_result(changed)
    changed = copy.deepcopy(result)
    changed["suite"] = "standard"
    validate_result(changed)
    for field in ("median_seconds", "rounds", "input_hash"):
        changed = copy.deepcopy(result)
        changed["cases"][0][field] = None
        with pytest.raises(ValueError, match=field):
            validate_result(changed)


def test_schema_version_and_states_after_d9():
    """The schema version stays 1. Only the cold state remains."""
    assert _api("benchmarks.fault_evaluation.run", "SCHEMA_VERSION") == 1
    assert _api("benchmarks.fault_evaluation.run", "STATES") == ("cold",)


def _row_of(case_id: str, **fields) -> dict:
    """Return a case row that has only the fields that cold_view reads."""
    model, operation, points, state = case_id.split(":")
    return {
        "case_id": case_id,
        "model": model,
        "operation": operation,
        "points": int(points),
        "state": state,
        **fields,
    }


def test_cold_view_turns_a_lone_warm_row_into_a_cold_row_and_drops_other_warm_rows():
    """A warm row with no cold row of the same case becomes cold. Other warm rows go (C-016)."""
    cold_view = _api("benchmarks.fault_evaluation.run", "cold_view")
    document = {
        "schema_version": 1,
        "cases": [
            _row_of("a:oldest_value:5000:cold", median_seconds=1.0),
            _row_of("a:oldest_value:5000:warm", median_seconds=9.0),
            _row_of("b:oldest_value:5000:warm", median_seconds=2.0),
            _row_of("c:model_ids:5000:cold", median_seconds=3.0),
        ],
        "provenance": {"python": "3.x"},
    }
    original = copy.deepcopy(document)
    view = cold_view(document)
    assert document == original, "cold_view must return a copy"
    assert view["provenance"] == document["provenance"]
    by_id = {row["case_id"]: row for row in view["cases"]}
    assert set(by_id) == {
        "a:oldest_value:5000:cold",
        "b:oldest_value:5000:cold",
        "c:model_ids:5000:cold",
    }
    assert by_id["a:oldest_value:5000:cold"]["median_seconds"] == 1.0
    assert by_id["b:oldest_value:5000:cold"]["median_seconds"] == 2.0
    assert {row["state"] for row in view["cases"]} == {"cold"}


def test_suite_run_case_makes_only_the_timed_call_and_never_reads_the_cache():
    """The suite case runs no priming call and no cache clear (D9). Warmup stays 0."""
    run_case = _api("benchmarks.fault_evaluation.bench_cases", "run_case")
    bench_module = importlib.import_module("benchmarks.fault_evaluation.bench_cases")
    events = []

    class FakeModel:
        @property
        def evaluation_cache(self):
            raise AssertionError("run_case must not read evaluation_cache")

    class TinyBenchmark:
        def __init__(self):
            self.extra_info = {}

        def pedantic(self, function, *, setup, rounds, warmup_rounds):
            self.warmup_rounds = warmup_rounds
            self.rounds = rounds
            for _ in range(rounds):
                setup()
                function()

    original = bench_module._OPERATIONS["model_ids"]
    bench_module._OPERATIONS["model_ids"] = lambda model, points: events.append("timed")
    try:
        benchmark = TinyBenchmark()
        case = {
            "case_id": "tiny:model_ids:2:cold",
            "model": "tiny",
            "state": "cold",
            "operation": "model_ids",
            "points": 2,
        }
        run_case(benchmark, FakeModel(), np.zeros((2, 3)), case, rounds=2)
    finally:
        bench_module._OPERATIONS["model_ids"] = original
    assert benchmark.warmup_rounds == 0
    assert benchmark.rounds == 2
    assert events == ["timed", "timed"]


def test_rss_conversion_normalizes_supported_native_units_to_bytes():
    """RSS unit conversion is explicit and rejects unknown platform units."""
    rss_to_bytes = _api("benchmarks.fault_evaluation.run", "rss_to_bytes")
    assert rss_to_bytes(2, "bytes") == 2
    assert rss_to_bytes(2, "kibibytes") == 2048
    with pytest.raises(ValueError, match="unit"):
        rss_to_bytes(2, "pages")
