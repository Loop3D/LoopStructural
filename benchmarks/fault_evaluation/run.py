"""Run the fault evaluation benchmarks and write one schema v1 result file.

Usage: ``python -m benchmarks.fault_evaluation.run --output PATH --rounds 1``

Time, counters and peak memory run in separate fresh processes. The counter
wrappers and the memory sampler then cannot change the timings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

from benchmarks.fault_evaluation.instrumentation import rss_to_bytes
from tests.fault_evaluation.models import model_names

SCHEMA_VERSION = 1
ROOT = Path(__file__).parents[2]
SUITE = Path(__file__).with_name("bench_cases.py")
PASS_KINDS = ("timing", "counters", "rss")
OPERATIONS = (
    "oldest_value",
    "all_group_values",
    "model_ids",
    "oldest_gradient",
    "model_gradient",
)
SIZES = (5_000, 50_000, 100_000, 500_000)
EXTRA_SIZE = 20_000
EXTRA_MODELS = ("plugin", "f1_u0")
STATES = ("cold", "warm")
_now = time.monotonic  # tests replace the clock

COUNTER_FIELDS = (
    "fault_calls_total",
    "fault_calls_top_level",
    "field_calls_total",
    "fault_points_total",
    "field_points_total",
)
RSS_FIELDS = ("peak_extra_rss_bytes", "rss_unit")
CASE_IDENTITY_FIELDS = ("case_id", "model", "operation", "points", "state")
CASE_FIELDS = (
    *CASE_IDENTITY_FIELDS,
    "rounds",
    "time_samples_seconds",
    "median_seconds",
    *COUNTER_FIELDS,
    *RSS_FIELDS,
    "code_revision",
    "dirty",
    "prerequisite_parent_revision",
    "prerequisite_fix_revision",
    "baseline_revision",
    "input_hash",
)
PROVENANCE_FIELDS = ("python", "numpy", "scipy", "loopstructural", "os", "cpu", "ram_bytes")

# The two revisions of the gradient prerequisite (plan D0).
PREREQUISITE_PARENT_REVISION = "b22517f984971b98f6d039f3dd3cde0b31b7600d"
PREREQUISITE_FIX_REVISION = "e9b72293de32344b7fe486bb24338c4efcd37cc8"
_DATA = ROOT / "tests" / "fault_evaluation" / "data"
_PLUGIN_DATA = ROOT / "specs" / "003-fault-evaluation-performance" / "data" / "plugin_model"


def enumerate_cases() -> list[dict[str, Any]]:
    """Return every benchmark case in a stable order."""
    sizes_by_model = {name: list(SIZES) for name in model_names()}
    for name in EXTRA_MODELS:
        sizes_by_model[name].append(EXTRA_SIZE)
    return [
        {
            "case_id": f"{name}:{operation}:{count}:{state}",
            "model": name,
            "operation": operation,
            "points": count,
            "state": state,
        }
        for name, sizes in sizes_by_model.items()
        for operation in OPERATIONS
        for count in sizes
        for state in STATES
    ]


# The 24 large cases of the standard suite (plan D2). They were drawn once with
# numpy.random.default_rng(2003). The pool for each size was the cases of f0_u0..f0_u3,
# f1_u0..f1_u3, stopped_fault, onlap_domain and faulted_fold (five operations, cold and
# warm), sorted by case id. For each size, in the order 50,000 then 100,000 and with one
# generator, rng.choice(len(pool), size=12, replace=False) chose the cases. The ids are
# sorted in the order of the full list. Do not draw them again: the text is the source.
STANDARD_LARGE_CASES: tuple[str, ...] = (
    "f0_u0:all_group_values:100000:warm",
    "f0_u0:oldest_gradient:100000:warm",
    "f0_u1:oldest_value:100000:warm",
    "f0_u1:model_ids:50000:warm",
    "f0_u1:model_gradient:100000:cold",
    "f0_u2:oldest_value:50000:warm",
    "f0_u3:model_ids:100000:warm",
    "f1_u0:all_group_values:100000:warm",
    "f1_u1:oldest_value:50000:cold",
    "f1_u1:model_gradient:100000:cold",
    "f1_u2:oldest_value:100000:warm",
    "f1_u2:oldest_gradient:50000:warm",
    "f1_u3:oldest_value:50000:cold",
    "f1_u3:all_group_values:100000:warm",
    "stopped_fault:model_ids:50000:cold",
    "stopped_fault:model_ids:50000:warm",
    "stopped_fault:oldest_gradient:50000:cold",
    "stopped_fault:oldest_gradient:50000:warm",
    "onlap_domain:model_ids:50000:cold",
    "onlap_domain:model_ids:50000:warm",
    "onlap_domain:oldest_gradient:100000:cold",
    "onlap_domain:model_gradient:100000:warm",
    "faulted_fold:oldest_value:50000:warm",
    "faulted_fold:all_group_values:100000:cold",
)
SUITE_NAMES = ("standard", "full")


def standard_cases() -> list[dict[str, Any]]:
    """Return the standard suite: all small cases and the 24 large cases, in full-list order."""
    large = set(STANDARD_LARGE_CASES)
    return [
        case
        for case in enumerate_cases()
        if case["points"] in (SIZES[0], EXTRA_SIZE) or case["case_id"] in large
    ]


def select_cases(suite: str) -> list[dict[str, Any]]:
    """Return the cases of the named suite. Raise ValueError for an unknown name."""
    if suite == "standard":
        return standard_cases()
    if suite == "full":
        return enumerate_cases()
    raise ValueError(f"Unknown suite {suite!r}. Use one of {SUITE_NAMES}.")


def validate_result(result: dict[str, Any]) -> None:
    """Raise ValueError if a result does not follow schema v1."""
    if result.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported schema_version {result.get('schema_version')!r}")
    provenance = result.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("Missing provenance")  # noqa: TRY004
    for field in PROVENANCE_FIELDS:
        if field not in provenance:
            raise ValueError(f"Missing provenance field {field}")
    if "pytest_benchmark" not in result:
        raise ValueError("Missing pytest_benchmark data")
    seen: set[str] = set()
    for case in result.get("cases", []):
        for field in CASE_FIELDS:
            if field not in case:
                raise ValueError(f"Case {case.get('case_id')!r} is missing field {field}")
        if case["case_id"] in seen:
            raise ValueError(f"Duplicate case_id {case['case_id']!r}")
        seen.add(case["case_id"])
        if case["rss_unit"] != "bytes":
            raise ValueError(
                f"Case {case['case_id']!r} has rss unit {case['rss_unit']!r}, not bytes"
            )
    if not seen:
        raise ValueError("The result has no case")


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _baseline_revision(code_revision: str) -> str:
    """Return the baseline revision of the golden manifest, if it exists."""
    try:
        manifest = json.loads((_DATA / "manifest.json").read_text(encoding="utf8"))
        return str(manifest["provenance"]["baseline_revision"])
    except (OSError, KeyError, ValueError):
        return code_revision


def _input_hash() -> str:
    """Hash the frozen inputs of all models."""
    digest = hashlib.sha256()
    paths = [_DATA / "recipes.json", _DATA / "fields.npz"]
    paths += sorted(path for path in _PLUGIN_DATA.glob("*") if path.is_file())
    for path in paths:
        digest.update(path.name.encode("utf8"))
        digest.update(path.read_bytes() if path.is_file() else b"missing")
    return digest.hexdigest()


def _ram_bytes() -> int:
    if sys.platform == "win32":
        import ctypes

        class _MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        status = _MemoryStatus()
        status.dwLength = ctypes.sizeof(_MemoryStatus)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
        return int(status.ullTotalPhys)
    try:
        return int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
    except (OSError, ValueError, AttributeError):
        return 0


def _provenance() -> dict[str, Any]:
    import numpy
    import scipy

    import LoopStructural

    return {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "loopstructural": getattr(LoopStructural, "__version__", "unknown"),
        "os": platform.platform(),
        "cpu": platform.processor() or platform.machine() or "unknown",
        "ram_bytes": _ram_bytes(),
    }


def _launch_pytest(suite: Path | str, args: list[str]) -> int:
    """Run pytest on the benchmark suite in a new Python process."""
    command = [sys.executable, "-m", "pytest", str(suite), *args]
    return subprocess.run(command, cwd=ROOT, check=False).returncode


# Paths that the runner writes itself. They never make the repository dirty.
_DIRTY_EXCLUDES: tuple[Path, ...] = ()


def _set_dirty_excludes(*paths: Path) -> None:
    """Set the paths that the dirty flag ignores. Call this before the first metadata call."""
    global _DIRTY_EXCLUDES
    _DIRTY_EXCLUDES = tuple(Path(path).resolve() for path in paths)
    cache_clear = getattr(_revision_metadata, "cache_clear", None)
    if cache_clear is not None:
        cache_clear()


def _is_dirty() -> bool:
    """Return True if the repository has changes, ignoring the paths of the runner."""
    root = ROOT.resolve()
    pathspecs = ["."]
    for path in _DIRTY_EXCLUDES:
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue  # The path is outside the repository.
        if relative.parts:
            pathspecs.append(":(exclude,literal)" + relative.as_posix())
    return bool(_git("status", "--porcelain", "--", *pathspecs))


@lru_cache(maxsize=1)
def _revision_metadata() -> dict[str, Any]:
    """Return the revision and input fields that every case row shares."""
    code_revision = _git("rev-parse", "HEAD")
    return {
        "code_revision": code_revision,
        "dirty": _is_dirty(),
        "prerequisite_parent_revision": PREREQUISITE_PARENT_REVISION,
        "prerequisite_fix_revision": PREREQUISITE_FIX_REVISION,
        "baseline_revision": _baseline_revision(code_revision),
        "input_hash": _input_hash(),
    }


def _case_row(case: dict[str, Any], entry: dict[str, Any], pass_kind: str) -> dict[str, Any]:
    stats = entry["stats"]
    info = entry["extra_info"]
    timing = pass_kind == "timing"
    row = {
        **{field: case[field] for field in CASE_IDENTITY_FIELDS},
        "rounds": int(stats["rounds"]),
        "time_samples_seconds": list(stats["data"]) if timing else [],
        "median_seconds": stats["median"] if timing else None,
        **dict.fromkeys(COUNTER_FIELDS, 0),
        "peak_extra_rss_bytes": 0,
        "rss_unit": "bytes",
        **_revision_metadata(),
    }
    if pass_kind == "counters":
        row.update({field: int(info["counters"][field]) for field in COUNTER_FIELDS})
    elif pass_kind == "rss":
        row["peak_extra_rss_bytes"] = rss_to_bytes(info["peak_extra_rss_bytes"], "bytes")
    return row


def _read_checkpoint(path: Path) -> dict[str, dict[str, Any]]:
    """Read the checkpoint lines by case id. Skip a line that a stop cut short."""
    entries: dict[str, dict[str, Any]] = {}
    if not path.is_file():
        return entries
    for line in path.read_text(encoding="utf8").splitlines():
        try:
            entry = json.loads(line)
            entries[entry["case_id"]] = entry
        except (ValueError, KeyError, TypeError):
            continue
    return entries


def _run_fresh_process_pass(
    pass_kind: str, cases: list[dict[str, Any]], *, output_path: Path, rounds: int
) -> dict[str, Any]:
    """Run one pass of the suite in a fresh process and write its result."""
    if pass_kind not in PASS_KINDS:
        raise ValueError(f"Unknown pass kind {pass_kind!r}")
    output_path = Path(output_path)
    cases_path = output_path.with_suffix(".cases.json")
    benchmark_path = output_path.with_suffix(".benchmark.json")
    checkpoint_path = output_path.with_suffix(".checkpoint.jsonl")
    # Cases that the checkpoint file already holds do not run again.
    entries = _read_checkpoint(checkpoint_path)
    pending = [case for case in cases if case["case_id"] not in entries]
    started = time.monotonic()
    print(
        f"[{pass_kind}] start: {len(pending)} of {len(cases)} cases to run",
        file=sys.stderr,
        flush=True,
    )
    if pending:
        cases_path.write_text(json.dumps(pending), encoding="utf8")
        benchmark_path.unlink(missing_ok=True)
        environment = {
            "FAULT_EVAL_PASS": pass_kind,
            "FAULT_EVAL_CASES": str(cases_path),
            "FAULT_EVAL_ROUNDS": str(rounds),
            "FAULT_EVAL_CHECKPOINT": str(checkpoint_path),
        }
        previous = {key: os.environ.get(key) for key in environment}
        os.environ.update(environment)
        try:
            args = ["--benchmark-json", str(benchmark_path), "--benchmark-save-data", "-v"]
            code = _launch_pytest(SUITE, args)
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
        if code != 0:
            raise RuntimeError(f"The {pass_kind} pass failed with exit code {code}")
        entries = _read_checkpoint(checkpoint_path)
        if benchmark_path.is_file():
            # A case without a checkpoint line comes from the raw pytest-benchmark file.
            raw = json.loads(benchmark_path.read_text(encoding="utf8"))
            for item in raw["benchmarks"]:
                entries.setdefault(item["extra_info"]["case_id"], item)
    missing = [case["case_id"] for case in cases if case["case_id"] not in entries]
    if missing:
        raise RuntimeError(f"The {pass_kind} pass has no result for {missing[:5]}")
    benchmarks = [
        {
            "extra_info": entries[case["case_id"]]["extra_info"],
            "stats": entries[case["case_id"]]["stats"],
        }
        for case in cases
    ]
    payload = {
        "schema_version": SCHEMA_VERSION,
        "pass": pass_kind,
        "cases": [_case_row(case, entries[case["case_id"]], pass_kind) for case in cases],
        "provenance": _provenance(),
        "pytest_benchmark": {"raw": {"benchmarks": benchmarks}} if pass_kind == "timing" else {},
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf8")
    print(
        f"[{pass_kind}] end: {len(cases)} cases, {time.monotonic() - started:.1f} s elapsed",
        file=sys.stderr,
        flush=True,
    )
    return payload


def _merge(
    passes: dict[str, dict[str, Any]], suite_wall_seconds: float | None = None
) -> dict[str, Any]:
    """Combine the timing, counter and memory passes by case id."""
    counters = {row["case_id"]: row for row in passes["counters"]["cases"]}
    memory = {row["case_id"]: row for row in passes["rss"]["cases"]}
    rows = []
    for row in passes["timing"]["cases"]:
        merged = dict(row)
        for field in COUNTER_FIELDS:
            merged[field] = counters[row["case_id"]][field]
        for field in RSS_FIELDS:
            merged[field] = memory[row["case_id"]][field]
        rows.append(merged)
    result = {
        "schema_version": SCHEMA_VERSION,
        "cases": rows,
        "provenance": passes["timing"]["provenance"],
        "pytest_benchmark": passes["timing"]["pytest_benchmark"],
    }
    if suite_wall_seconds is not None:
        result["suite_wall_seconds"] = suite_wall_seconds
    return result


def _prepare_work_dir(
    parser: argparse.ArgumentParser, work_dir: Path, *, resume: bool, rounds: int
) -> dict[str, Any]:
    """Check the work directory and return the meta data of the run.

    Stop with a clear error if a resume does not fit the saved run, or if a new
    run would overwrite saved data.
    """
    meta_path = work_dir / "meta.json"
    current = _revision_metadata()
    if resume and not meta_path.is_file() and work_dir.exists() and any(work_dir.iterdir()):
        parser.error(
            f"Cannot resume: the work directory {work_dir} holds data but no meta.json. "
            "Delete the directory to start again."
        )
    if not resume:
        if work_dir.exists() and any(work_dir.iterdir()):
            parser.error(
                f"The work directory {work_dir} already holds data. "
                "Pass --resume to continue that run, or delete the directory."
            )
    elif meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf8"))
        for key, value in (
            ("code_revision", current["code_revision"]),
            ("dirty", current["dirty"]),
            ("rounds", rounds),
        ):
            if meta.get(key) != value:
                parser.error(
                    f"Cannot resume: saved {key} {meta.get(key)!r} differs from "
                    f"current {value!r}. Delete {work_dir} to start again."
                )
        return meta
    work_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "code_revision": current["code_revision"],
        "dirty": current["dirty"],
        "rounds": rounds,
        "wall_seconds": dict.fromkeys(PASS_KINDS, 0.0),
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf8")
    return meta


def main(argv: list[str] | None = None) -> int:
    """Run all three passes and write the merged result."""
    parser = argparse.ArgumentParser(description="Run the fault evaluation benchmarks.")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Directory for the checkpoint files (default: the --output path plus .work).",
    )
    parser.add_argument(
        "--resume", action="store_true", help="Continue a stopped run in the work directory."
    )
    parser.add_argument(
        "--suite",
        choices=SUITE_NAMES,
        default="standard",
        help="Cases to run: the standard suite (default) or the full suite.",
    )
    options = parser.parse_args(argv)
    if options.rounds < 1:
        parser.error("--rounds must be at least 1")
    work_dir = options.work_dir or options.output.with_name(options.output.name + ".work")
    # The work directory and the output file are not changes of the code.
    _set_dirty_excludes(work_dir, options.output)
    meta = _prepare_work_dir(parser, work_dir, resume=options.resume, rounds=options.rounds)
    cases = select_cases(options.suite)
    passes = {}
    for pass_kind in PASS_KINDS:
        output_path = work_dir / f"{pass_kind}.json"
        if options.resume and output_path.is_file():
            passes[pass_kind] = json.loads(output_path.read_text(encoding="utf8"))
            print(f"[{pass_kind}] already complete, skipped", file=sys.stderr, flush=True)
            continue
        started = _now()
        try:
            passes[pass_kind] = _run_fresh_process_pass(
                pass_kind, cases, output_path=output_path, rounds=options.rounds
            )
        finally:
            # A failed or stopped pass also adds its time.
            meta["wall_seconds"][pass_kind] += _now() - started
            (work_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf8")
    result = _merge(passes, suite_wall_seconds=float(sum(meta["wall_seconds"].values())))
    validate_result(result)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
