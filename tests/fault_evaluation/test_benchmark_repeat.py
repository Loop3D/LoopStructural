"""Tests for the repeat of short cases in the timing pass of the benchmark runner (SC-007).

All tests use a fake launcher. No real pytest process and no real model runs.
"""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path

import pytest

run_module = importlib.import_module("benchmarks.fault_evaluation.run")

CHECKPOINT_ENV = "FAULT_EVAL_CHECKPOINT"
# Median of each case in seconds. Cases 0 and 1 are short. Cases 2 and 3 are not.
MEDIANS = {0: 0.2, 1: 0.9, 2: 1.0, 3: 3.0}
SHORT_IDS = ["tiny:0", "tiny:1"]


def _cases() -> list[dict]:
    return [
        {
            "case_id": f"tiny:{index}",
            "model": "tiny",
            "operation": "oldest_value",
            "points": 2 + index,
            "state": "cold",
        }
        for index in MEDIANS
    ]


class _Launcher:
    """A fake pytest launcher. It writes checkpoint lines and records each launch."""

    def __init__(self, clock: dict | None = None):
        self.launches: list[dict] = []
        self.clock = clock
        # Launch number (1 based) and number of cases to write before the stop.
        self.stop_at: tuple[int, int] | None = None

    def __call__(self, suite, args):
        cases = json.loads(Path(os.environ["FAULT_EVAL_CASES"]).read_text(encoding="utf8"))
        rounds = int(os.environ["FAULT_EVAL_ROUNDS"])
        self.launches.append(
            {
                "pass": os.environ["FAULT_EVAL_PASS"],
                "cases": [case["case_id"] for case in cases],
                "rounds": rounds,
            }
        )
        if self.clock is not None:
            self.clock["now"] += 10.0
        limit = None
        if self.stop_at and self.stop_at[0] == len(self.launches):
            limit = self.stop_at[1]
        for index, case in enumerate(cases):
            if limit is not None and index >= limit:
                return 1
            number = int(case["case_id"].split(":")[1])
            median = MEDIANS[number]
            # The samples show how many rounds produced them.
            line = {
                "case_id": case["case_id"],
                "stats": {
                    "rounds": rounds,
                    "data": [median + 0.001 * step for step in range(rounds)],
                    "median": median,
                },
                "extra_info": {
                    **case,
                    "counters": dict.fromkeys(run_module.COUNTER_FIELDS, number + 1),
                    "peak_extra_rss_bytes": 1024 * (number + 1),
                },
                "elapsed_seconds": 0.5,
            }
            with open(os.environ[CHECKPOINT_ENV], "a", encoding="utf8") as handle:
                handle.write(json.dumps(line) + "\n")
        return 0


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


def _prepare_main(monkeypatch, clock: dict | None = None):
    monkeypatch.setattr(run_module, "enumerate_cases", _cases)
    clock = clock if clock is not None else {"now": 0.0}
    monkeypatch.setattr(run_module, "_now", lambda: clock["now"])


def _argv(tmp_path, *extra):
    return [
        "--output",
        str(tmp_path / "result.json"),
        "--work-dir",
        str(tmp_path / "work"),
        "--suite",
        "full",  # the fake enumerate_cases list is the full suite
        *extra,
    ]


def _rows(result: dict) -> dict:
    return {row["case_id"]: row for row in result["cases"]}


def test_constants_are_defined():
    """The runner names the limit of a short case and the number of repeat rounds."""
    assert run_module.SHORT_CASE_SECONDS == 1.0
    assert run_module.SHORT_CASE_ROUNDS == 5


def test_short_cases_are_repeated_with_five_rounds(monkeypatch, tmp_path):
    """The second step runs only the cases with a median below 1.0 s, with 5 rounds."""
    launcher = _Launcher()
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    payload = run_module._run_fresh_process_pass(
        "timing", _cases(), output_path=tmp_path / "timing.json", rounds=1
    )
    assert launcher.launches == [
        {"pass": "timing", "cases": [f"tiny:{i}" for i in MEDIANS], "rounds": 1},
        {"pass": "timing", "cases": SHORT_IDS, "rounds": run_module.SHORT_CASE_ROUNDS},
    ]
    rows = _rows(payload)
    for case_id in SHORT_IDS:
        assert rows[case_id]["rounds"] == 5
        assert len(rows[case_id]["time_samples_seconds"]) == 5
    # A case of exactly 1.0 s or more keeps its one round.
    for case_id in ("tiny:2", "tiny:3"):
        assert rows[case_id]["rounds"] == 1
        assert len(rows[case_id]["time_samples_seconds"]) == 1
    # The new samples replace the old ones also in the raw benchmark data.
    raw = payload["pytest_benchmark"]["raw"]["benchmarks"]
    by_id = {item["extra_info"]["case_id"]: item for item in raw}
    assert len(by_id["tiny:0"]["stats"]["data"]) == 5
    assert len(by_id["tiny:2"]["stats"]["data"]) == 1
    # The checkpoint file keeps the later line for a case.
    checkpoint = run_module._read_checkpoint(tmp_path / "timing.checkpoint.jsonl")
    assert checkpoint["tiny:0"]["stats"]["rounds"] == 5


def test_no_repeat_when_first_step_has_enough_rounds(monkeypatch, tmp_path):
    """A first step with 5 rounds or more needs no second step."""
    launcher = _Launcher()
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    payload = run_module._run_fresh_process_pass(
        "timing", _cases(), output_path=tmp_path / "timing.json", rounds=5
    )
    assert len(launcher.launches) == 1
    assert {row["rounds"] for row in payload["cases"]} == {5}


@pytest.mark.parametrize("pass_kind", ["counters", "rss"])
def test_counters_and_memory_passes_are_not_repeated(monkeypatch, tmp_path, pass_kind):
    """Only the timing pass has a second step."""
    launcher = _Launcher()
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    run_module._run_fresh_process_pass(
        pass_kind, _cases(), output_path=tmp_path / f"{pass_kind}.json", rounds=1
    )
    assert [launch["rounds"] for launch in launcher.launches] == [1]


def test_main_repeats_only_in_timing_pass(monkeypatch, tmp_path):
    """The merged result holds 5 rounds for short cases and one launch per other pass."""
    _prepare_main(monkeypatch)
    launcher = _Launcher()
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    assert run_module.main(_argv(tmp_path)) == 0
    assert [(launch["pass"], launch["rounds"]) for launch in launcher.launches] == [
        ("timing", 1),
        ("timing", 5),
        ("counters", 1),
        ("rss", 1),
    ]
    assert launcher.launches[2]["cases"] == [f"tiny:{i}" for i in MEDIANS]
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf8"))
    run_module.validate_result(result)
    rows = _rows(result)
    assert [rows[case_id]["rounds"] for case_id in rows] == [5, 5, 1, 1]
    assert rows["tiny:0"]["median_seconds"] == pytest.approx(0.2)
    # The counters and memory values still come from their own passes.
    assert rows["tiny:0"][run_module.COUNTER_FIELDS[0]] == 1
    assert rows["tiny:1"]["peak_extra_rss_bytes"] == 2048


def test_stop_between_the_steps_does_not_lose_or_repeat_twice(monkeypatch, tmp_path):
    """A stop before the second step, then a resume: each short case repeats once."""
    _prepare_main(monkeypatch)
    launcher = _Launcher()
    launcher.stop_at = (2, 0)  # the second launch stops before it writes a line
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    with pytest.raises(RuntimeError, match="exit code"):
        run_module.main(_argv(tmp_path))
    assert not (tmp_path / "work" / "timing.json").exists()

    assert run_module.main(_argv(tmp_path, "--resume")) == 0
    timing_launches = [x for x in launcher.launches if x["pass"] == "timing"]
    # The first step is complete, so the resume runs only the second step, once.
    assert [launch["cases"] for launch in timing_launches[1:]] == [SHORT_IDS, SHORT_IDS]
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf8"))
    rows = _rows(result)
    assert [rows[case_id]["rounds"] for case_id in rows] == [5, 5, 1, 1]


def test_stop_inside_the_second_step_repeats_only_missing_cases(monkeypatch, tmp_path):
    """A case with 5 rounds in the checkpoint is not repeated again after a resume."""
    _prepare_main(monkeypatch)
    launcher = _Launcher()
    launcher.stop_at = (2, 1)  # the second launch writes tiny:0 and then stops
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    with pytest.raises(RuntimeError, match="exit code"):
        run_module.main(_argv(tmp_path))

    assert run_module.main(_argv(tmp_path, "--resume")) == 0
    timing_launches = [x for x in launcher.launches if x["pass"] == "timing"]
    assert [launch["cases"] for launch in timing_launches] == [
        [f"tiny:{i}" for i in MEDIANS],
        SHORT_IDS,
        ["tiny:1"],
    ]
    assert timing_launches[2]["rounds"] == run_module.SHORT_CASE_ROUNDS
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf8"))
    assert [row["rounds"] for row in result["cases"]] == [5, 5, 1, 1]


def test_resume_of_complete_timing_pass_does_not_repeat(monkeypatch, tmp_path):
    """A resume after the timing pass finished launches nothing for timing."""
    _prepare_main(monkeypatch)
    launcher = _Launcher()
    launcher.stop_at = (3, 0)  # the counters pass stops
    monkeypatch.setattr(run_module, "_launch_pytest", launcher)
    with pytest.raises(RuntimeError):
        run_module.main(_argv(tmp_path))
    assert run_module.main(_argv(tmp_path, "--resume")) == 0
    assert [x["pass"] for x in launcher.launches].count("timing") == 2


def test_repeat_adds_its_time_to_suite_wall_seconds(monkeypatch, tmp_path):
    """Each launch takes 10 s on the fake clock: 2 timing launches, counters and memory."""
    clock = {"now": 0.0}
    _prepare_main(monkeypatch, clock)
    monkeypatch.setattr(run_module, "_launch_pytest", _Launcher(clock))
    assert run_module.main(_argv(tmp_path)) == 0
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf8"))
    assert result["suite_wall_seconds"] == pytest.approx(40.0)
