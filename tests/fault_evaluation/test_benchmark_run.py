"""Contract tests for opt-in fault evaluation benchmark tooling."""

from __future__ import annotations

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
    """All required models/operations/sizes have cold and warm stable cases."""
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
        for state in ("cold", "warm")
    }
    observed = {(item["model"], item["operation"], item["points"], item["state"]) for item in cases}
    assert required <= observed
    assert {item["model"] for item in cases} == models

    for operation in operations:
        for state in ("cold", "warm"):
            assert ("plugin", operation, 20_000, state) in observed
            assert ("f1_u0", operation, 20_000, state) in observed


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


def test_runner_separates_passes_and_writes_merged_case_metrics(monkeypatch, tmp_path):
    """Three mocked passes merge into schema output with timing raw statistics."""
    run_main = _api("benchmarks.fault_evaluation.run", "main")
    run_module = importlib.import_module(run_main.__module__)
    cases = [
        {
            "case_id": "tiny:cold",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2,
            "state": "cold",
        },
        {
            "case_id": "tiny:warm",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2,
            "state": "warm",
        },
    ]
    monkeypatch.setattr(run_module, "enumerate_cases", lambda: cases)
    calls = []

    def fake_pass(pass_kind, selected_cases, *, output_path, rounds):
        calls.append((pass_kind, [item["case_id"] for item in selected_cases], rounds))
        payload = _tiny_result(pass_kind, selected_cases)
        Path(output_path).write_text(json.dumps(payload), encoding="utf8")
        return payload

    monkeypatch.setattr(run_module, "_run_fresh_process_pass", fake_pass)
    monkeypatch.setattr(run_module, "_launch_pytest", lambda *args: 0)
    output = tmp_path / "merged.json"
    assert run_main(["--output", str(output), "--rounds", "2", "--suite", "full"]) == 0
    assert [call[0] for call in calls] == ["timing", "counters", "rss"]
    assert all(call[1] == ["tiny:cold", "tiny:warm"] and call[2] == 2 for call in calls)
    merged = json.loads(output.read_text(encoding="utf8"))
    assert merged["schema_version"] == 1
    assert [row["case_id"] for row in merged["cases"]] == ["tiny:cold", "tiny:warm"]
    assert all(row["fault_calls_total"] == 2 for row in merged["cases"])
    assert all(row["field_points_total"] == 7 for row in merged["cases"])
    assert all(row["peak_extra_rss_bytes"] == 4096 for row in merged["cases"])
    assert merged["pytest_benchmark"]["raw"]["stats"]["mean"] == 0.015


@pytest.mark.parametrize("state,expected_hits,expected_misses", [("cold", 0, 2), ("warm", 2, 0)])
def test_suite_run_case_controls_warmup_and_cold_cache(state, expected_hits, expected_misses):
    """The suite case owns pedantic warmup and untimed cache preparation."""
    run_case = _api("benchmarks.fault_evaluation.bench_cases", "run_case")
    benchmark = None

    class FakeModel:
        def __init__(self):
            self.cached = False
            self.hits = 0
            self.misses = 0
            self.warm_primes = 0
            self.evaluation_cache = FakeCache(self)

        def evaluate_model(self, points):
            return self._evaluate(points)

        def _evaluate(self, points):
            if not benchmark.in_timing:
                assert state == "warm", "cold targets must not receive untimed warmup calls"
                assert self.warm_primes < 2, "warm setup should prime once before each timed round"
                self.warm_primes += 1
            else:
                if self.cached:
                    self.hits += 1
                else:
                    self.misses += 1
            self.cached = True
            return np.zeros(len(points))

    class FakeCache:
        def __init__(self, model):
            self.model = model

        def clear(self):
            assert not benchmark.in_timing
            self.model.cached = False

    model = FakeModel()
    points = np.zeros((2, 3))

    class TinyBenchmark:
        def __init__(self):
            self.in_timing = False
            self.extra_info = {}

        def pedantic(self, function, *, setup, rounds, warmup_rounds):
            benchmark.in_timing = False
            self.warmup_rounds = warmup_rounds
            self.rounds = rounds
            for _ in range(rounds):
                setup()
                self.in_timing = True
                try:
                    function()
                finally:
                    self.in_timing = False

    benchmark = TinyBenchmark()
    case = {
        "case_id": f"tiny:{state}",
        "model": "tiny",
        "state": state,
        "operation": "model_ids",
        "points": 2,
    }
    run_case(benchmark, model, points, case, rounds=2)
    assert benchmark.warmup_rounds == 0
    assert benchmark.rounds == 2
    assert model.warm_primes == (2 if state == "warm" else 0)
    assert model.hits == expected_hits
    assert model.misses == expected_misses


def test_rss_conversion_normalizes_supported_native_units_to_bytes():
    """RSS unit conversion is explicit and rejects unknown platform units."""
    rss_to_bytes = _api("benchmarks.fault_evaluation.run", "rss_to_bytes")
    assert rss_to_bytes(2, "bytes") == 2
    assert rss_to_bytes(2, "kibibytes") == 2048
    with pytest.raises(ValueError, match="unit"):
        rss_to_bytes(2, "pages")
