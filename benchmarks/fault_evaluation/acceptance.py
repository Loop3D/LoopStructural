"""Check the numerical performance targets SC-002 to SC-010 on schema v1 results.

Usage::

    python -m benchmarks.fault_evaluation.acceptance --baseline PATH --candidate PATH \
        --repeat PATH --output PATH [--memory-points 10000000]

``check_acceptance`` compares a baseline result with a candidate result. It
reports one row for each target. A failed target fails the whole result. It is
never a warning. The candidate holds three more inputs: ``repeat`` (a second
full result, for SC-007), ``suite_wall_seconds`` (for SC-008) and
``gradient_decision`` (the OP1 rule).

``--memory-points`` is a local mode for the large-memory check (FR-010). It runs
each cache mode in a fresh subprocess. Do not run it in pytest or in CI.
"""

# ruff: noqa: ANN401  (the wrapped objects are untyped model and fixture objects)

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
GRADIENT_RULE = "op1-one-restoration-per-fault-set-and-point-batch"
SHORT_CASE_SECONDS = 1.0
NOISE_LIMIT = 1.10
SUITE_WALL_LIMIT_SECONDS = 900.0
MEMORY_MODES = ("cold", "warm", "off")
# FR-010: 3 floats of 8 bytes for each point and each fault set.
RESTORED_BYTES_PER_POINT = 3 * 8
_TOLERANCE = 1e-9

OLD_F4U3 = "f4_u3:oldest_value:5000:cold"
ALL_F4U3 = "f4_u3:all_group_values:5000:cold"
ALL_F4U0 = "f4_u0:all_group_values:5000:cold"
WARM_F4U3 = "f4_u3:all_group_values:5000:warm"
# SC-010 uses evaluate_model, which is the benchmark operation model_ids.
IDS_COLD_F4U3 = "f4_u3:model_ids:5000:cold"
IDS_WARM_F4U3 = "f4_u3:model_ids:5000:warm"
PLUGIN_OLD = "plugin:oldest_value:20000:cold"
SINGLE_F1U0 = "f1_u0:oldest_value:20000:cold"
_PROVENANCE_FIELDS = ("python", "numpy", "scipy", "os", "cpu", "ram_bytes")
_CASE_PROVENANCE_FIELDS = (
    "input_hash",
    "baseline_revision",
    "prerequisite_parent_revision",
    "prerequisite_fix_revision",
)


def _row(passed: bool, value: Any, limit: Any, reason: str | None) -> dict[str, Any]:
    """Return one criterion row. A passed row has no reason."""
    return {"passed": passed, "value": value, "limit": limit, "reason": None if passed else reason}


def _index(document: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """Return the cases by id (first one wins) and the ids that occur twice."""
    cases: dict[str, dict[str, Any]] = {}
    duplicates: list[str] = []
    for case in document.get("cases", []):
        case_id = case["case_id"]
        if case_id in cases:
            duplicates.append(case_id)
        else:
            cases[case_id] = case
    return cases, duplicates


def _check_cases(
    base: dict[str, dict[str, Any]],
    new: dict[str, dict[str, Any]],
    base_duplicates: list[str],
    new_duplicates: list[str],
) -> dict[str, Any]:
    """Check that both results hold the same cases, each one once."""
    problems = [f"Duplicate case {case_id!r} in the baseline" for case_id in base_duplicates]
    problems += [f"Duplicate case {case_id!r} in the candidate" for case_id in new_duplicates]
    problems += [
        f"Case {case_id!r} is missing from the candidate" for case_id in base if case_id not in new
    ]
    problems += [
        f"Case {case_id!r} is missing from the baseline" for case_id in new if case_id not in base
    ]
    return {"passed": not problems, "problems": problems}


def _check_provenance(
    base: dict[str, Any],
    new: dict[str, Any],
    base_cases: dict[str, dict[str, Any]],
    new_cases: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Check schema version, machine fields and the input identity of each case."""
    problems = []
    for label, document in (("baseline", base), ("candidate", new)):
        if document.get("schema_version") != SCHEMA_VERSION:
            problems.append(
                f"The {label} has schema_version {document.get('schema_version')!r}, "
                f"expected {SCHEMA_VERSION}"
            )
    base_provenance = base.get("provenance", {})
    new_provenance = new.get("provenance", {})
    for field in _PROVENANCE_FIELDS:
        if base_provenance.get(field) != new_provenance.get(field):
            problems.append(
                f"Provenance field {field} differs: baseline {base_provenance.get(field)!r}, "
                f"candidate {new_provenance.get(field)!r}"
            )
    for case_id, base_case in base_cases.items():
        new_case = new_cases.get(case_id)
        if new_case is None:
            continue
        for field in _CASE_PROVENANCE_FIELDS:
            if base_case.get(field) != new_case.get(field):
                problems.append(
                    f"Field {field} differs for case {case_id!r}: baseline "
                    f"{base_case.get(field)!r}, candidate {new_case.get(field)!r}"
                )
    return {"passed": not problems, "problems": problems}


def _missing(cases: dict[str, dict[str, Any]], *case_ids: str) -> str | None:
    """Return a reason that names the first case id that is absent."""
    for case_id in case_ids:
        if case_id not in cases:
            return f"Case {case_id!r} is missing"
    return None


def _time_ratio_row(
    base: dict[str, dict[str, Any]],
    new: dict[str, dict[str, Any]],
    case_id: str,
    limit: float,
    calls: int,
) -> dict[str, Any]:
    """Check time ratio against the baseline and the exact count of fault calls."""
    reason = _missing(base, case_id) or _missing(new, case_id)
    if reason:
        return _row(False, None, limit, reason)
    base_median = base[case_id]["median_seconds"]
    if base_median <= 0:
        return _row(False, None, limit, f"Case {case_id!r} has a baseline median of zero")
    ratio = new[case_id]["median_seconds"] / base_median
    problems = []
    if ratio > limit + _TOLERANCE:
        problems.append(f"time ratio {ratio:.4f} of case {case_id!r} is above {limit}")
    if new[case_id]["fault_calls_total"] != calls:
        problems.append(
            f"case {case_id!r} makes {new[case_id]['fault_calls_total']} fault calls, not {calls}"
        )
    return _row(not problems, ratio, limit, "; ".join(problems))


def _sc004(new: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """SC-004: time with 3 unconformities over time with 0, both from the candidate."""
    limit = 1.5
    reason = _missing(new, ALL_F4U3, ALL_F4U0)
    if reason:
        return _row(False, None, limit, reason)
    plain = new[ALL_F4U0]["median_seconds"]
    if plain <= 0:
        return _row(False, None, limit, f"Case {ALL_F4U0!r} has a median of zero")
    ratio = new[ALL_F4U3]["median_seconds"] / plain
    return _row(ratio <= limit + _TOLERANCE, ratio, limit, f"ratio {ratio:.4f} is above {limit}")


def _sc006(base: dict[str, dict[str, Any]], new: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """SC-006: the peak extra memory of every case is at most 120 % of the baseline."""
    limit = 1.20
    worst = 0.0
    problems = []
    for case_id, base_case in base.items():
        if case_id not in new:
            problems.append(f"case {case_id!r} is missing from the candidate")
            continue
        base_rss = base_case["peak_extra_rss_bytes"]
        new_rss = new[case_id]["peak_extra_rss_bytes"]
        if base_rss <= 0:
            if new_rss > 0:
                problems.append(f"case {case_id!r} has a baseline memory of zero")
            continue
        ratio = new_rss / base_rss
        worst = max(worst, ratio)
        if ratio > limit + _TOLERANCE:
            problems.append(f"case {case_id!r} memory ratio {ratio:.4f} is above {limit}")
    return _row(not problems, worst, limit, "; ".join(problems))


def _check_noise(
    base: dict[str, dict[str, Any]],
    new: dict[str, dict[str, Any]],
    repeat: dict[str, Any] | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """SC-007 (plan D2b). Return the criterion row and the list of short cases."""
    if repeat is None:
        return _row(False, None, NOISE_LIMIT, "The repeat run is missing"), []
    repeat_cases, repeat_duplicates = _index(repeat)
    problems = [f"duplicate case {case_id!r} in the repeat run" for case_id in repeat_duplicates]
    short_cases = []
    worst = 0.0
    for case_id, base_case in base.items():
        if case_id not in new or case_id not in repeat_cases:
            where = "candidate" if case_id not in new else "repeat run"
            problems.append(f"case {case_id!r} is missing from the {where}")
            continue
        medians = (
            base_case["median_seconds"],
            new[case_id]["median_seconds"],
            repeat_cases[case_id]["median_seconds"],
        )
        if min(medians) <= 0:
            problems.append(f"case {case_id!r} has a median of zero")
            continue
        noise = max(medians[1] / medians[2], medians[2] / medians[1])
        if medians[0] < SHORT_CASE_SECONDS:
            short_cases.append({"case_id": case_id, "noise_ratio": noise})
            continue
        worst = max(worst, noise)
        if noise >= NOISE_LIMIT:
            problems.append(f"case {case_id!r} has noise ratio {noise:.4f}")
    return _row(not problems, worst, NOISE_LIMIT, "; ".join(problems)), short_cases


def _sc008(new: dict[str, Any]) -> dict[str, Any]:
    """SC-008: the standard suite runs in less than 900 s."""
    wall = new.get("suite_wall_seconds")
    if wall is None:
        return _row(False, None, SUITE_WALL_LIMIT_SECONDS, "suite_wall_seconds is missing")
    return _row(
        wall < SUITE_WALL_LIMIT_SECONDS,
        wall,
        SUITE_WALL_LIMIT_SECONDS,
        f"suite time {wall} s is not below {SUITE_WALL_LIMIT_SECONDS} s",
    )


def _sc009(base: dict[str, dict[str, Any]], new: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """SC-009: the time for each point of fault application is at most 50 % of the baseline."""
    limit = 0.50
    reason = _missing(base, SINGLE_F1U0) or _missing(new, SINGLE_F1U0)
    if reason:
        return _row(False, None, limit, reason)
    per_point = []
    for cases in (base, new):
        points = cases[SINGLE_F1U0]["fault_points_total"]
        if points <= 0:
            return _row(False, None, limit, f"Case {SINGLE_F1U0!r} has no fault points")
        per_point.append(cases[SINGLE_F1U0]["median_seconds"] / points)
    if per_point[0] <= 0:
        return _row(False, None, limit, f"Case {SINGLE_F1U0!r} has a baseline time of zero")
    ratio = per_point[1] / per_point[0]
    return _row(ratio <= limit + _TOLERANCE, ratio, limit, f"ratio {ratio:.4f} is above {limit}")


def _sc010(new: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """SC-010: a second evaluate_model call makes 0 fault calls and takes at most 10 %.

    evaluate_model is the benchmark operation ``model_ids``. The cold case of the
    same document is the reference. SC-003 keeps the operation all_group_values.
    """
    limit = 0.10
    reason = _missing(new, IDS_WARM_F4U3, IDS_COLD_F4U3)
    if reason:
        return _row(False, None, limit, reason)
    cold = new[IDS_COLD_F4U3]["median_seconds"]
    if cold <= 0:
        return _row(False, None, limit, f"Case {IDS_COLD_F4U3!r} has a median of zero")
    ratio = new[IDS_WARM_F4U3]["median_seconds"] / cold
    problems = []
    if ratio > limit + _TOLERANCE:
        problems.append(f"warm ratio {ratio:.4f} of case {IDS_WARM_F4U3!r} is above {limit}")
    calls = new[IDS_WARM_F4U3]["fault_calls_total"]
    if calls != 0:
        problems.append(f"case {IDS_WARM_F4U3!r} makes {calls} fault calls")
    return _row(not problems, ratio, limit, "; ".join(problems))


def check_acceptance(base: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Check SC-002 to SC-010 for a baseline and a candidate result.

    The candidate holds ``repeat``, ``suite_wall_seconds`` and ``gradient_decision``.
    The result is passed only if the case check, the provenance check, the
    gradient decision and every criterion pass.
    """
    base_cases, base_duplicates = _index(base)
    new_cases, new_duplicates = _index(new)
    cases = _check_cases(base_cases, new_cases, base_duplicates, new_duplicates)
    provenance = _check_provenance(base, new, base_cases, new_cases)
    decision = new.get("gradient_decision")
    gradient = {
        "passed": decision == GRADIENT_RULE,
        "expected": GRADIENT_RULE,
        "value": decision,
    }
    noise, short_cases = _check_noise(base_cases, new_cases, new.get("repeat"))
    criteria = {
        "SC-002": _time_ratio_row(base_cases, new_cases, OLD_F4U3, 0.40, 4),
        "SC-003": _time_ratio_row(base_cases, new_cases, ALL_F4U3, 0.15, 4),
        "SC-004": _sc004(new_cases),
        "SC-005": _time_ratio_row(base_cases, new_cases, PLUGIN_OLD, 0.50, 6),
        "SC-006": _sc006(base_cases, new_cases),
        "SC-007": noise,
        "SC-008": _sc008(new),
        "SC-009": _sc009(base_cases, new_cases),
        "SC-010": _sc010(new_cases),
    }
    passed = (
        cases["passed"]
        and provenance["passed"]
        and gradient["passed"]
        and all(row["passed"] for row in criteria.values())
    )
    return {
        "passed": passed,
        "criteria": criteria,
        "cases": cases,
        "provenance": provenance,
        "gradient_decision": gradient,
        "short_cases": short_cases,
    }


# --- Large-memory mode (local only) -------------------------------------------------------


def _cache_fields(model: Any) -> dict[str, Any]:
    """Return the cache fields of a memory record.

    The cache interface (plan D4) has the read-only property ``retained_bytes``:
    the exact sum of ``nbytes`` of the stored arrays. This function reads only
    that property. If the model has a cache but the property is missing or is not
    an int, the record gets null bytes and an ``error`` text. Then the memory
    mode fails. It never records null without an error.
    """
    cache = getattr(model, "evaluation_cache", None)
    if cache is None:
        return {"cache_available": False, "retained_cache_bytes": None}
    value = getattr(cache, "retained_bytes", None)
    if isinstance(value, int) and not isinstance(value, bool):
        return {"cache_available": True, "retained_cache_bytes": value}
    return {
        "cache_available": True,
        "retained_cache_bytes": None,
        "error": f"The evaluation cache has no int property retained_bytes (got {value!r})",
    }


def _distinct_fault_sets(model: Any) -> int:
    """Return the number of distinct ordered fault sets over the features of the model.

    A fault set is the tuple of fault names of one feature. Features with no
    faults do not count. A model with no faults gives 0.
    """
    sets = set()
    for feature in getattr(model, "features", None) or []:
        names = tuple(fault.name for fault in getattr(feature, "faults", None) or [])
        if names:
            sets.add(names)
    return len(sets)


def _fr010_accounting(
    points: int, fault_sets: int, batches: dict[int, int], peak: int
) -> dict[str, Any]:
    """Return the FR-010 accounting fields of a memory record.

    ``allowed_bytes`` is 24 bytes x points x distinct fault sets.
    ``restored_bytes_observed`` is 24 bytes x the sum of (batch size x count) over
    the batches that the FaultSegment.apply_to_points wrapper counted.
    ``peak_over_allowed`` is the peak extra RSS over ``allowed_bytes``. It is null
    when ``allowed_bytes`` is 0.

    This is a report only. The exit code does not depend on it. The accepted OP1
    rule measures each actual restored point and excludes unchanged baseline work
    arrays. RSS cannot separate these two parts.
    """
    allowed = RESTORED_BYTES_PER_POINT * points * fault_sets
    restored_points = sum(size * count for size, count in batches.items())
    return {
        "distinct_fault_sets": fault_sets,
        "allowed_bytes": allowed,
        "restored_bytes_observed": RESTORED_BYTES_PER_POINT * restored_points,
        "peak_over_allowed": peak / allowed if allowed else None,
    }


def _memory_worker(model_name: str, points: int, mode: str) -> dict[str, Any]:
    """Run one large evaluation in this process and return its memory record."""
    from benchmarks.fault_evaluation.bench_cases import _points
    from benchmarks.fault_evaluation.instrumentation import RssSampler
    from LoopStructural.modelling.features.fault import FaultSegment
    from tests.fault_evaluation.models import build_model

    model = build_model(model_name)
    cache = getattr(model, "evaluation_cache", None)
    fault_sets = _distinct_fault_sets(model)
    coordinates = _points(model, points)
    if cache is not None:
        cache.clear()
        cache.enabled = mode != "off"
    if mode == "warm" and cache is not None:
        model.evaluate_model(coordinates)
    batches: Counter = Counter()
    original = FaultSegment.apply_to_points

    def apply_to_points(fault: Any, locations: Any, *args: Any, **kwargs: Any) -> Any:
        batches[len(locations)] += 1
        return original(fault, locations, *args, **kwargs)

    FaultSegment.apply_to_points = apply_to_points
    sampler = RssSampler()
    sampler.start()
    try:
        model.evaluate_model(coordinates)
    finally:
        peak = sampler.stop()
        FaultSegment.apply_to_points = original
    return {
        "model": model_name,
        "mode": mode,
        "points": points,
        **_cache_fields(model),
        "peak_extra_rss_bytes": peak,
        "restored_batch_sizes": {str(size): count for size, count in sorted(batches.items())},
        **_fr010_accounting(points, fault_sets, dict(batches), peak),
        "input_bytes": coordinates.nbytes,
        "fr010": {
            "bytes_per_point_and_fault_set": RESTORED_BYTES_PER_POINT,
            "bytes_per_fault_set": RESTORED_BYTES_PER_POINT * points,
        },
    }


def run_memory_mode(points: int, models: list[str]) -> dict[str, Any]:
    """Run each model and cache mode in a fresh subprocess. Return the records."""
    records = []
    for model_name in models:
        for mode in MEMORY_MODES:
            command = [
                sys.executable,
                "-m",
                "benchmarks.fault_evaluation.acceptance",
                "--memory-worker",
                model_name,
                str(points),
                mode,
            ]
            done = subprocess.run(command, capture_output=True, text=True, check=False)
            if done.returncode != 0:
                records.append(
                    {"model": model_name, "mode": mode, "points": points, "error": done.stderr}
                )
                continue
            records.append(json.loads(done.stdout.strip().splitlines()[-1]))
    return {"points": points, "records": records}


def main(argv: list[str] | None = None) -> int:
    """Run the acceptance check, the large-memory mode, or both. Return 0 on success."""
    parser = argparse.ArgumentParser(description="Check the fault evaluation performance targets.")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--repeat", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--memory-points",
        type=int,
        help="Run the large-memory check with this point count (documented value: 10000000).",
    )
    parser.add_argument(
        "--memory-model",
        action="append",
        help="Frozen model for the memory check. Can repeat. Default: f1_u0.",
    )
    parser.add_argument("--memory-worker", nargs=3, metavar=("MODEL", "POINTS", "MODE"))
    options = parser.parse_args(argv)
    if options.memory_worker:
        model_name, points, mode = options.memory_worker
        print(json.dumps(_memory_worker(model_name, int(points), mode)))
        return 0
    if options.output is None:
        parser.error("--output is required")
    paths = (options.baseline, options.candidate, options.repeat)
    check = all(path is not None for path in paths)
    if not check and options.memory_points is None:
        parser.error("--baseline, --candidate and --repeat are required")
    if not check and any(path is not None for path in paths):
        parser.error("--baseline, --candidate and --repeat must be given together")
    result: dict[str, Any] = {}
    passed = True
    if check:
        base = json.loads(options.baseline.read_text(encoding="utf8"))
        new = json.loads(options.candidate.read_text(encoding="utf8"))
        new["repeat"] = json.loads(options.repeat.read_text(encoding="utf8"))
        result = check_acceptance(base, new)
        passed = result["passed"]
        for name, row in result["criteria"].items():
            status = "PASS" if row["passed"] else "FAIL"
            print(
                f"{name} {status} value={row['value']} limit={row['limit']} {row['reason'] or ''}"
            )
        print(f"overall {'PASS' if passed else 'FAIL'}")
    if options.memory_points is not None:
        result["memory"] = run_memory_mode(options.memory_points, options.memory_model or ["f1_u0"])
        # Only a record with an error fails the run. Do not gate the exit code on RSS:
        # the OP1 rule counts each actual restored point and excludes unchanged
        # baseline work arrays, and RSS cannot separate these two parts.
        failed = [record for record in result["memory"]["records"] if "error" in record]
        passed = passed and not failed
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
