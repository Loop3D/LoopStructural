"""Tests for progress, checkpoint and resume of the fault evaluation benchmark runner.

All tests use a fake launcher. No real pytest process and no real model runs.
"""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

run_module = importlib.import_module("benchmarks.fault_evaluation.run")
bench_module = importlib.import_module("benchmarks.fault_evaluation.bench_cases")

CHECKPOINT_ENV = "FAULT_EVAL_CHECKPOINT"


def _cases(count: int = 4) -> list[dict]:
    return [
        {
            "case_id": f"tiny:{index}",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2 + index,
            "state": "cold",
        }
        for index in range(count)
    ]


def _fake_launcher(stop_after: int | None = None, launches: list | None = None):
    """Return a fake pytest launcher that writes checkpoint lines like the suite."""

    def launch(suite, args):
        cases = json.loads(Path(os.environ["FAULT_EVAL_CASES"]).read_text(encoding="utf8"))
        if launches is not None:
            launches.append([case["case_id"] for case in cases])
        checkpoint = os.environ.get(CHECKPOINT_ENV)
        for index, case in enumerate(cases):
            if stop_after is not None and index >= stop_after:
                return 1
            number = int(case["case_id"].split(":")[1])
            line = {
                "case_id": case["case_id"],
                "stats": {
                    "rounds": 2,
                    "data": [0.01 * (number + 1), 0.02 * (number + 1)],
                    "median": 0.015 * (number + 1),
                },
                "extra_info": {
                    **case,
                    "counters": dict.fromkeys(run_module.COUNTER_FIELDS, number + 1),
                    "peak_extra_rss_bytes": 1024 * (number + 1),
                },
                "elapsed_seconds": 0.5,
            }
            with open(checkpoint, "a", encoding="utf8") as handle:
                handle.write(json.dumps(line) + "\n")
        return 0

    return launch


@pytest.fixture(autouse=True)
def _fixed_environment(monkeypatch):
    metadata = {
        "code_revision": "rev1",
        "dirty": False,
        "prerequisite_parent_revision": "parent",
        "prerequisite_fix_revision": "fix",
        "baseline_revision": "baseline",
        "input_hash": "0" * 64,
    }
    monkeypatch.setattr(run_module, "_revision_metadata", lambda: metadata)
    monkeypatch.setattr(
        run_module,
        "_provenance",
        lambda: {
            "python": "3.x",
            "numpy": "2.x",
            "scipy": "1.x",
            "loopstructural": "2.x",
            "os": "test",
            "cpu": "test cpu",
            "ram_bytes": 1,
        },
    )
    return metadata


@pytest.mark.parametrize("pass_kind", ["timing", "counters", "rss"])
def test_stopped_and_resumed_pass_equals_straight_pass(monkeypatch, tmp_path, pass_kind):
    """A pass that stops half way and resumes gives the same rows as a straight pass."""
    cases = _cases()
    monkeypatch.setattr(run_module, "_launch_pytest", _fake_launcher())
    (tmp_path / "straight").mkdir()
    straight = run_module._run_fresh_process_pass(
        pass_kind, cases, output_path=tmp_path / "straight" / f"{pass_kind}.json", rounds=2
    )

    launches: list = []
    output = tmp_path / "resumed" / f"{pass_kind}.json"
    output.parent.mkdir()
    monkeypatch.setattr(run_module, "_launch_pytest", _fake_launcher(2, launches))
    with pytest.raises(RuntimeError, match="exit code"):
        run_module._run_fresh_process_pass(pass_kind, cases, output_path=output, rounds=2)
    monkeypatch.setattr(run_module, "_launch_pytest", _fake_launcher(None, launches))
    resumed = run_module._run_fresh_process_pass(pass_kind, cases, output_path=output, rounds=2)

    assert launches == [["tiny:0", "tiny:1", "tiny:2", "tiny:3"], ["tiny:2", "tiny:3"]]
    assert resumed["cases"] == straight["cases"]
    assert resumed["pytest_benchmark"].keys() == straight["pytest_benchmark"].keys()


def test_pass_uses_verbose_progress_and_sets_checkpoint_env(monkeypatch, tmp_path, capsys):
    """The pass shows each case (-v), names the checkpoint file and prints start and end."""
    seen = {}

    def launch(suite, args):
        seen["args"] = args
        seen["checkpoint"] = os.environ.get(CHECKPOINT_ENV)
        return _fake_launcher()(suite, args)

    monkeypatch.setattr(run_module, "_launch_pytest", launch)
    run_module._run_fresh_process_pass(
        "counters", _cases(2), output_path=tmp_path / "counters.json", rounds=1
    )
    assert "-v" in seen["args"] and "-q" not in seen["args"]
    assert seen["checkpoint"] and seen["checkpoint"].startswith(str(tmp_path))
    assert CHECKPOINT_ENV not in os.environ
    lines = [line for line in capsys.readouterr().err.splitlines() if "counters" in line]
    assert len(lines) == 2


def _prepare_main(monkeypatch, cases, clock_step=10.0):
    monkeypatch.setattr(run_module, "enumerate_cases", lambda: cases)
    ticks = iter(range(10_000))
    monkeypatch.setattr(run_module, "_now", lambda: next(ticks) * clock_step)


def _argv(tmp_path, *extra):
    return [
        "--output",
        str(tmp_path / "result.json"),
        "--work-dir",
        str(tmp_path / "work"),
        *extra,
    ]


def test_resume_is_refused_for_other_revision_or_dirty_flag(
    monkeypatch, tmp_path, _fixed_environment
):
    """--resume stops with a clear error when the code revision or dirty flag changed."""
    _prepare_main(monkeypatch, _cases(2))
    monkeypatch.setattr(run_module, "_launch_pytest", _fake_launcher())
    assert run_module.main(_argv(tmp_path)) == 0

    _fixed_environment["code_revision"] = "rev2"
    with pytest.raises(SystemExit) as error:
        run_module.main(_argv(tmp_path, "--resume"))
    assert error.value.code != 0
    _fixed_environment["code_revision"] = "rev1"
    _fixed_environment["dirty"] = True
    with pytest.raises(SystemExit) as error:
        run_module.main(_argv(tmp_path, "--resume"))
    assert error.value.code != 0
    _fixed_environment["dirty"] = False
    assert run_module.main(_argv(tmp_path, "--resume")) == 0


def test_non_empty_work_dir_without_resume_is_refused(monkeypatch, tmp_path, capsys):
    """Without --resume, a work directory that holds data is never overwritten."""
    _prepare_main(monkeypatch, _cases(2))
    monkeypatch.setattr(run_module, "_launch_pytest", _fake_launcher())
    work = tmp_path / "work"
    work.mkdir()
    assert run_module.main(_argv(tmp_path)) == 0  # an empty directory is accepted
    before = sorted(path.name for path in work.iterdir())
    with pytest.raises(SystemExit) as error:
        run_module.main(_argv(tmp_path))
    assert error.value.code != 0
    assert "--resume" in capsys.readouterr().err
    assert sorted(path.name for path in work.iterdir()) == before


def test_suite_wall_seconds_accumulates_across_resume(monkeypatch, tmp_path):
    """The result holds the sum of the wall seconds of all passes, also across resumes."""
    _prepare_main(monkeypatch, _cases(2), clock_step=10.0)
    state = {"fail": True}
    real_launch = _fake_launcher()

    def launch(suite, args):
        if state["fail"] and os.environ["FAULT_EVAL_PASS"] == "counters":
            return 1
        return real_launch(suite, args)

    monkeypatch.setattr(run_module, "_launch_pytest", launch)
    with pytest.raises(RuntimeError):
        run_module.main(_argv(tmp_path))  # timing done (10 s), counters fails (10 s)
    assert not (tmp_path / "result.json").exists()
    state["fail"] = False
    assert run_module.main(_argv(tmp_path, "--resume")) == 0  # counters + rss: 20 s
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf8"))
    assert result["suite_wall_seconds"] == pytest.approx(40.0)
    run_module.validate_result(result)
    del result["suite_wall_seconds"]
    run_module.validate_result(result)


def test_no_checkpoint_is_written_without_the_environment_variable(monkeypatch, tmp_path):
    """The suite hook writes nothing when FAULT_EVAL_CHECKPOINT is absent."""
    monkeypatch.delenv(CHECKPOINT_ENV, raising=False)
    monkeypatch.chdir(tmp_path)
    stats = SimpleNamespace(rounds=2, data=[0.1, 0.2], median=0.15)
    benchmark = SimpleNamespace(stats=SimpleNamespace(stats=stats), extra_info={"case_id": "a"})
    bench_module.write_checkpoint(benchmark, {"case_id": "a"}, 1.5)
    assert list(tmp_path.iterdir()) == []

    target = tmp_path / "sub.jsonl"
    monkeypatch.setenv(CHECKPOINT_ENV, str(target))
    bench_module.write_checkpoint(benchmark, {"case_id": "a"}, 1.5)
    line = json.loads(target.read_text(encoding="utf8"))
    assert line["case_id"] == "a"
    assert line["stats"] == {"rounds": 2, "data": [0.1, 0.2], "median": 0.15}
    assert line["extra_info"]["case_id"] == "a"
    assert line["elapsed_seconds"] == 1.5
