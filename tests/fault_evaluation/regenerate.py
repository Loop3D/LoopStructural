"""Validate or explicitly regenerate frozen fault-evaluation references."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import scipy
import LoopStructural

from tests.fault_evaluation.golden import capture
from tests.fault_evaluation.models import build_model, model_names, sample_points

REFERENCE_DIR = Path(__file__).parent / "data"
_PROVENANCE = {
    "prerequisite_parent_revision": "b22517f984971b98f6d039f3dd3cde0b31b7600d",
    "prerequisite_fix_revision": "e9b72293de32344b7fe486bb24338c4efcd37cc8",
}


def _revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).parents[2], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_paths() -> dict[str, Path]:
    root = Path(__file__).parents[2]
    paths = {
        "golden.py": Path(__file__).with_name("golden.py"),
        "regenerate.py": Path(__file__),
        "models.py": Path(__file__).with_name("models.py"),
        "recipes.json": REFERENCE_DIR / "recipes.json",
        "fields.npz": REFERENCE_DIR / "fields.npz",
        "geological_feature.py": root / "LoopStructural/modelling/features/_geological_feature.py",
        "base_geological_feature.py": root
        / "LoopStructural/modelling/features/_base_geological_feature.py",
        "geological_model.py": root / "LoopStructural/modelling/core/geological_model.py",
        "structured_grid.py": root
        / "packages/loop_common/src/loop_common/supports/_3d_structured_grid.py",
    }
    for source_root in (
        root / "LoopStructural",
        root / "packages/loop_common/src",
        root / "packages/loop_interpolation/src",
    ):
        if source_root.is_dir():
            for path in source_root.rglob("*.py"):
                paths[str(path.relative_to(root)).replace("\\", "/")] = path
    return paths


def _source_hashes() -> dict[str, str]:
    return {
        name: _hash(path) if path.is_file() else "missing" for name, path in _source_paths().items()
    }


def _git_status() -> str:
    try:
        return subprocess.check_output(
            ["git", "status", "--short"], cwd=Path(__file__).parents[2], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _plugin_hashes() -> dict[str, str]:
    plugin_dir = (
        Path(__file__).parents[2] / "specs/003-fault-evaluation-performance/data/plugin_model"
    )
    return {path.name: _hash(path) for path in sorted(plugin_dir.glob("*")) if path.is_file()}


def _load(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def _box_size(model: Any) -> float:
    box = model.bounding_box
    return float(np.max(np.asarray(box.maximum) - np.asarray(box.origin)))


def _manifest() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "provenance": {
            **_PROVENANCE,
            "baseline_revision": _revision(),
            "baseline_dirty": bool(_git_status()),
            "dirty_details": _git_status(),
            "field_generation_revision": "0a7930e2",
            "source_hashes": _source_hashes(),
            "input_hashes": _plugin_hashes(),
        },
        "dependencies": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "LoopStructural": getattr(LoopStructural, "__version__", "unknown"),
            "platform": platform.platform(),
        },
        "seeds": {"interior": 1, "count": 128},
        "point_recipe": {
            "candidate_seed_offset": 1701,
            "candidate_count": 2048,
            "feature_offset_fraction": 1e-4,
            "categories": [
                "interior",
                "box_faces",
                "fault_surfaces",
                "unit_boundaries",
                "unconformities",
            ],
        },
        "models": {},
    }


def _validate() -> int:
    manifest_path = REFERENCE_DIR / "manifest.json"
    if not manifest_path.is_file():
        return 1
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf8"))
        for name in model_names():
            path = REFERENCE_DIR / f"{name}.npz"
            entry = manifest["models"][name]
            if not path.is_file() or _hash(path) != entry["sha256"]:
                return 1
            model = build_model(name)
            expected = _load(path)
            actual = capture(model, expected["points"])
            assert_matches(name, actual, expected, _box_size(model))
    except Exception:
        return 1
    return 0


def assert_matches(
    name: str, actual: dict[str, np.ndarray], expected: dict[str, np.ndarray], box_size: float
) -> None:
    """Imported lazily here to keep the public module API convenient."""
    from tests.fault_evaluation.golden import assert_matches as compare

    compare(name, actual, expected, box_size)


def main(argv: list[str] | None = None) -> int:
    """Validate references by default; write files only with explicit --write."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="write reference archives and manifest"
    )
    args = parser.parse_args(argv)
    if not args.write:
        return _validate()

    manifest_path = REFERENCE_DIR / "manifest.json"
    if manifest_path.is_file():
        old = json.loads(manifest_path.read_text(encoding="utf8"))
        old_hashes = old.get("provenance", {}).get("source_hashes")
        old_inputs = old.get("provenance", {}).get("input_hashes")
        if (old_hashes and old_hashes != _source_hashes()) or (
            old_inputs and old_inputs != _plugin_hashes()
        ):
            raise RuntimeError(
                "Source hashes differ from the frozen baseline; refusing to replace references"
            )
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    manifest = _manifest()
    for name in model_names():
        print(f"Capturing golden outputs: {name}", flush=True)
        model = build_model(name)
        points, point_categories = _reference_points(model, name, 128, seed=1)
        points_path = REFERENCE_DIR / f"{name}.npz"
        outputs = capture(model, points)
        np.savez_compressed(points_path, **outputs)
        manifest["models"][name] = {
            "sha256": _hash(points_path),
            "point_count": int(len(points)),
            "point_categories": point_categories,
            "box_size": _box_size(model),
            "transform": np.asarray(
                getattr(model.bounding_box, "_world_to_local", np.eye(4))
            ).tolist(),
            "surface_levels": sorted(
                {
                    key.split("@", 1)[1].split("::", 1)[0]
                    for key in outputs
                    if key.startswith("surface::")
                }
            ),
            "summary": {key: np.asarray(value).shape for key, value in outputs.items()},
        }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf8")
    return 0


def _reference_points(
    model: Any, name: str, count: int, seed: int
) -> tuple[np.ndarray, dict[str, int]]:
    """Add deterministic near-feature and near-boundary samples to interior points."""
    interior = np.asarray(sample_points(name, count, seed=seed), dtype=float)
    box = model.bounding_box
    origin = np.asarray(box.origin, dtype=float)
    maximum = np.asarray(box.maximum, dtype=float)
    size = float(np.max(maximum - origin))
    offset = max(size * 1e-4, np.finfo(float).eps)
    faces = []
    for axis in range(3):
        for side, coordinate in ((0, origin[axis] + offset), (1, maximum[axis] - offset)):
            point = origin + (maximum - origin) * np.array([0.413, 0.537, 0.619])
            point[axis] = coordinate
            faces.append(point)
    category_points: dict[str, list[np.ndarray]] = {
        "interior": list(interior),
        "box_faces": faces,
        "fault_surfaces": [],
        "unit_boundaries": [],
        "unconformities": [],
    }
    rng = np.random.default_rng(seed + 1701)
    candidates = rng.random((2048, 3)) * (maximum - origin) + origin

    def near_level(feature: Any, level: float) -> list[np.ndarray]:
        values = np.asarray(feature.evaluate_value(candidates), dtype=float).reshape(-1)
        finite = np.isfinite(values)
        if not np.any(finite):
            return []
        indices = np.flatnonzero(finite)
        chosen = indices[np.argsort(np.abs(values[finite] - level))[:4]]
        output = []
        for index in chosen:
            point = candidates[index].copy()
            for _ in range(10):
                value = float(np.asarray(feature.evaluate_value(point[None, :])).reshape(-1)[0])
                if not np.isfinite(value):
                    break
                gradient = _value_gradient(feature, point, size)
                norm_sq = float(np.dot(gradient, gradient))
                if not np.isfinite(norm_sq) or norm_sq <= 0:
                    break
                step = ((value - level) / norm_sq) * gradient
                step_norm = np.linalg.norm(step)
                if step_norm > size * 0.25:
                    step *= size * 0.25 / step_norm
                point = np.clip(point - step, origin, maximum)
                if abs(value - level) <= max(1e-8 * max(abs(level), 1.0), 1e-8):
                    break
            final_value = float(np.asarray(feature.evaluate_value(point[None, :])).reshape(-1)[0])
            if not np.isfinite(final_value) or abs(final_value - level) > max(
                1e-5 * max(abs(level), 1.0), 1e-5
            ):
                continue
            gradient = _value_gradient(feature, point, size)
            norm = np.linalg.norm(gradient)
            if not np.isfinite(norm) or norm == 0:
                continue
            direction = gradient / norm
            output.extend(
                (
                    np.clip(point + direction * offset, origin, maximum),
                    np.clip(point - direction * offset, origin, maximum),
                )
            )
        return output

    column = getattr(model, "stratigraphic_column", None)
    groups = (
        list(column.get_groups()) if column is not None and hasattr(column, "get_groups") else []
    )
    for group in groups:
        feature = next(
            (item for item in model.features if getattr(item, "name", None) == group.name), None
        )
        if feature is None:
            continue
        for unit in group.units:
            for level in (unit.min(), unit.max()):
                if np.isfinite(level):
                    found = near_level(feature, float(level))
                    category_points["unit_boundaries"].extend(found)
    for unconformity in model.features:
        value = getattr(unconformity, "value", None)
        parent = getattr(unconformity, "parent", None)
        if value is not None and parent is not None:
            category_points["unconformities"].extend(near_level(unconformity, float(value)))
    for fault in getattr(model, "faults", []):
        found = near_level(fault[0], 0.0)
        category_points["fault_surfaces"].extend(found)
    combined = np.vstack(
        [np.asarray(items, dtype=float).reshape(-1, 3) for items in category_points.values()]
    )
    rounded = np.round(combined, decimals=10)
    _, indices = np.unique(rounded, axis=0, return_index=True)
    selected = np.sort(indices)
    selected_keys = {tuple(row) for row in rounded[selected]}
    categories = {
        category: len({tuple(np.round(point, decimals=10)) for point in items} & selected_keys)
        for category, items in category_points.items()
    }
    return combined[selected], categories


def _value_gradient(feature: Any, point: np.ndarray, box_size: float) -> np.ndarray:
    """Estimate a scalar field direction using values only, without legacy refinement."""
    step = max(box_size * 1e-5, np.finfo(float).eps * 100)
    samples = np.repeat(np.asarray(point, dtype=float)[None, :], 6, axis=0)
    for axis in range(3):
        samples[2 * axis, axis] += step
        samples[2 * axis + 1, axis] -= step
    values = np.asarray(feature.evaluate_value(samples), dtype=float).reshape(-1)
    return np.asarray([(values[2 * axis] - values[2 * axis + 1]) / (2 * step) for axis in range(3)])


if __name__ == "__main__":
    sys.exit(main())
