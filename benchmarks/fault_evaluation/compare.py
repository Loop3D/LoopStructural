"""Compare two schema v1 benchmark results.

Usage: ``python -m benchmarks.fault_evaluation.compare BASE NEW --output PATH``

A time ratio is the new median time divided by the baseline median time. A
ratio above the noise limit is marked as slower. A case that does not match
between the two results is an error. The tool never skips a case.
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
NOISE_LIMIT = 1.10
COUNT_FIELDS = (
    "fault_calls_total",
    "field_calls_total",
    "fault_points_total",
    "field_points_total",
)
_MODEL_PATTERN = re.compile(r"^f(\d+)_u(\d+)$")


def _ratio(base: float, new: float) -> tuple[float | None, str]:
    """Return the new/base ratio and its status. A zero baseline has no ratio."""
    if base == 0:
        return None, "undefined_zero_baseline"
    return new / base, "ok"


def _check_schema(document: dict[str, Any], label: str) -> None:
    """Raise ValueError if a result does not use schema v1."""
    if document.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"The {label} result has unsupported schema_version "
            f"{document.get('schema_version')!r}"
        )


def _index_cases(document: dict[str, Any], label: str) -> dict[str, dict[str, Any]]:
    """Return the cases of one result by case id. Reject a duplicate id."""
    indexed: dict[str, dict[str, Any]] = {}
    for case in document["cases"]:
        case_id = case["case_id"]
        if case_id in indexed:
            raise ValueError(f"The {label} result has a duplicate case_id {case_id!r}")
        indexed[case_id] = case
    return indexed


def _case_row(base: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Return the comparison row for two matching cases."""
    ratio, status = _ratio(base["median_seconds"], new["median_seconds"])
    counts = {}
    for field in COUNT_FIELDS:
        count_ratio, count_status = _ratio(base[field], new[field])
        counts[field] = {
            "base": base[field],
            "new": new[field],
            "ratio": count_ratio,
            "ratio_status": count_status,
        }
    return {
        "case_id": base["case_id"],
        "model": base["model"],
        "operation": base["operation"],
        "points": base["points"],
        "state": base["state"],
        "base_median_seconds": base["median_seconds"],
        "new_median_seconds": new["median_seconds"],
        "ratio": ratio,
        "ratio_status": status,
        "slower": ratio is not None and ratio > NOISE_LIMIT,
        "counts": counts,
    }


def _scaling_table(groups: dict[int, list[float]], key: str) -> list[dict[str, Any]]:
    """Return one row per group with the median of its valid ratios."""
    return [
        {
            key: count,
            "ratio": statistics.median(ratios) if ratios else None,
            "cases": len(ratios),
        }
        for count, ratios in sorted(groups.items())
    ]


def _scaling(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Group the time ratios by fault count (F) and unconformity count (U).

    Only models named like ``f{F}_u{U}`` take part. The fault table keeps
    U = 0. The unconformity table keeps F = 0.
    """
    faults: dict[int, list[float]] = {}
    units: dict[int, list[float]] = {}
    for row in rows:
        match = _MODEL_PATTERN.match(row["model"])
        if match is None:
            continue
        fault_count, unit_count = int(match.group(1)), int(match.group(2))
        if unit_count == 0:
            faults.setdefault(fault_count, [])
            if row["ratio"] is not None:
                faults[fault_count].append(row["ratio"])
        if fault_count == 0:
            units.setdefault(unit_count, [])
            if row["ratio"] is not None:
                units[unit_count].append(row["ratio"])
    return {
        "fault": _scaling_table(faults, "fault_count"),
        "unconformity": _scaling_table(units, "unconformity_count"),
    }


def _provenance_warnings(base: dict[str, Any], new: dict[str, Any]) -> list[str]:
    """Return one warning for each provenance field that differs."""
    base_provenance = base["provenance"]
    new_provenance = new["provenance"]
    warnings = []
    for field in sorted(set(base_provenance) | set(new_provenance)):
        if base_provenance.get(field) != new_provenance.get(field):
            warnings.append(
                f"Provenance field {field} differs: base {base_provenance.get(field)!r}, "
                f"new {new_provenance.get(field)!r}"
            )
    return warnings


def compare_results(base: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Compare two schema v1 results case by case.

    Raises
    ------
    ValueError
        If a case is missing, duplicated, has a different point count, or has
        a different cold/warm state in the two results.
    """
    _check_schema(base, "base")
    _check_schema(new, "new")
    base_cases = _index_cases(base, "base")
    new_cases = _index_cases(new, "new")
    for case_id in base_cases:
        if case_id not in new_cases:
            raise ValueError(f"Case {case_id!r} is missing from the new result")
    for case_id in new_cases:
        if case_id not in base_cases:
            raise ValueError(f"Case {case_id!r} is missing from the base result")
    rows = []
    for case_id, base_case in base_cases.items():
        new_case = new_cases[case_id]
        if base_case["points"] != new_case["points"]:
            raise ValueError(
                f"Case {case_id!r} has different points: "
                f"base {base_case['points']!r}, new {new_case['points']!r}"
            )
        if base_case["state"] != new_case["state"]:
            raise ValueError(
                f"Case {case_id!r} has a different state: "
                f"base {base_case['state']!r}, new {new_case['state']!r}"
            )
        rows.append(_case_row(base_case, new_case))
    return {
        "schema_version": SCHEMA_VERSION,
        "noise_limit": NOISE_LIMIT,
        "cases": rows,
        "scaling": _scaling(rows),
        "warnings": _provenance_warnings(base, new),
    }


def main(argv: list[str] | None = None) -> int:
    """Compare two result files and write the comparison as JSON."""
    parser = argparse.ArgumentParser(description="Compare two fault evaluation benchmark results.")
    parser.add_argument("base", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    options = parser.parse_args(argv)
    base = json.loads(options.base.read_text(encoding="utf8"))
    new = json.loads(options.new.read_text(encoding="utf8"))
    try:
        result = compare_results(base, new)
    except ValueError as error:
        parser.error(str(error))
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
