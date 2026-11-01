"""Capture and compare deterministic fault-evaluation outputs."""

from __future__ import annotations

from typing import Any

import numpy as np


def _coordinates(model: Any) -> list[Any]:
    result: list[Any] = []
    seen: set[int] = set()
    pending = list(model.features)
    while pending:
        feature = pending.pop(0)
        if id(feature) in seen:
            continue
        seen.add(id(feature))
        result.append(feature)
        pending.extend(list(getattr(feature, "features", None) or []))
    return result


def _faulted_features(model: Any) -> list[tuple[Any, list[Any]]]:
    """Return one representative feature for every distinct ordered fault set."""
    result: list[tuple[Any, list[Any]]] = []
    keys: set[tuple[int, ...]] = set()
    for feature in _coordinates(model):
        faults = list(getattr(feature, "faults", None) or [])
        key = tuple(id(fault) for fault in faults)
        if faults and key not in keys:
            keys.add(key)
            result.append((feature, faults))
    return result


def _mesh_area(vertices: np.ndarray, triangles: np.ndarray) -> float:
    if not len(triangles):
        return 0.0
    tri = vertices[triangles]
    return float(
        np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum() / 2
    )


def _groups(model: Any) -> list[Any]:
    column = getattr(model, "stratigraphic_column", None)
    return list(column.get_groups()) if column is not None and hasattr(column, "get_groups") else []


def capture(model: Any, points: np.ndarray) -> dict[str, np.ndarray]:
    """Capture flat numeric model, feature, restoration and surface outputs."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    result: dict[str, np.ndarray] = {"points": points.copy()}
    ids = model.evaluate_model(points.copy())
    result["model_ids"] = np.asarray(ids)
    coordinates = _coordinates(model)
    for feature in coordinates:
        if not callable(getattr(feature, "evaluate_value", None)):
            continue
        name = str(getattr(feature, "name", type(feature).__name__))
        result[f"feature::{name}::value"] = np.asarray(feature.evaluate_value(points.copy()))
        gradient = getattr(feature, "evaluate_gradient", None)
        if callable(gradient):
            result[f"feature::{name}::gradient"] = np.asarray(gradient(points.copy()))
    for coordinate, faults in _faulted_features(model):
        names = ">".join(str(getattr(fault, "name", type(fault).__name__)) for fault in faults)
        restored = points.copy()
        if getattr(coordinate, "faults_enabled", True):
            for fault in faults:
                restored = fault.apply_to_points(restored)
        result[f"restored::{names}"] = np.asarray(restored)

    bounds = getattr(model, "bounding_box", None)
    from copy import deepcopy

    for group in _groups(model):
        # Stratigraphic columns can retain legacy group declarations whose
        # feature is absent from the loaded model (notably plugin columns).
        # Only capture surfaces for features that actually belong to the model.
        feature = next((f for f in model.features if getattr(f, "name", None) == group.name), None)
        if feature is None or not callable(getattr(feature, "surfaces", None)):
            continue
        group_levels = sorted(
            {float(v) for unit in group.units for v in (unit.min(), unit.max()) if np.isfinite(v)}
        )
        if not group_levels:
            continue
        box = deepcopy(bounds)
        box.nsteps = np.array([20, 20, 20])
        for level in group_levels:

            def record_expected_empty() -> None:
                """Accept emptiness only when the sampled field proves it."""
                grid_points = box.regular_grid(local=False, order="F")
                samples = np.asarray(feature.evaluate_value(grid_points))
                finite = samples[np.isfinite(samples)]
                mask_fn = getattr(feature, "_calculate_mask", None)
                visible = mask_fn(grid_points) if callable(mask_fn) else None
                fully_eroded = visible is not None and not np.any(visible)
                outside_range = finite.size and (level < finite.min() or level > finite.max())
                if not fully_eroded and not outside_range:
                    raise ValueError(f"Empty surface for {group.name} at {level:g}")
                _store_surface(result, feature.name, level, None, None)

            meshes = feature.surfaces([level], box)
            if len(meshes) > 1:
                raise ValueError(
                    f"Unexpected surface count for {group.name} at {level:g}: {len(meshes)}"
                )
            mesh = meshes[0] if meshes else None
            if mesh is None:
                record_expected_empty()
                continue
            vertices = np.asarray(getattr(mesh, "vertices", np.empty((0, 3))))
            triangles = np.asarray(getattr(mesh, "triangles", np.empty((0, 3), dtype=int)))
            if not len(vertices):
                record_expected_empty()
                continue
            _store_surface(result, feature.name, level, vertices, triangles)

    distances = np.empty((len(points), 0), dtype=float)
    for group in _groups(model):
        feature = next((f for f in model.features if getattr(f, "name", None) == group.name), None)
        key = f"feature::{group.name}::value"
        values = result.get(key)
        if feature is None or values is None:
            continue
        for unit in group.units:
            for boundary in (unit.min(), unit.max()):
                if not np.isfinite(boundary):
                    continue
                column_values = np.abs(np.asarray(values).reshape(-1) - boundary)
                column_values[~np.isfinite(column_values)] = np.inf
                distances = np.column_stack((distances, column_values))
    result["unit_boundary_distances"] = distances
    return result


def _store_surface(
    result: dict[str, np.ndarray],
    name: str,
    level: float,
    vertices: np.ndarray | None,
    triangles: np.ndarray | None,
) -> None:
    """Store a surface summary, including a valid empty erosion result."""
    key = f"surface::{name}@{level:g}"
    count = 0 if vertices is None else len(vertices)
    result[f"{key}::vertex_count"] = np.asarray([count], dtype=np.int64)
    result[f"{key}::bounds"] = (
        np.vstack((vertices.min(axis=0), vertices.max(axis=0)))
        if count
        else np.zeros((2, 3), dtype=float)
    )
    area = 0.0 if vertices is None or triangles is None else _mesh_area(vertices, triangles)
    result[f"{key}::area"] = np.asarray([area])


def assert_matches(
    name: str, actual: dict[str, np.ndarray], expected: dict[str, np.ndarray], box_size: float
) -> None:
    """Compare captures with the D1 tolerances and useful point diagnostics."""
    if set(actual) != set(expected):
        raise AssertionError(
            f"{name} key mismatch: actual-only={sorted(set(actual)-set(expected))}, expected-only={sorted(set(expected)-set(actual))}"
        )
    atol = 1e-8 * float(box_size)
    for key in expected:
        left, right = np.asarray(actual[key]), np.asarray(expected[key])
        if left.shape != right.shape:
            raise AssertionError(f"{name} {key} shape mismatch: {left.shape} != {right.shape}")
        if key == "model_ids":
            equal = left == right
            if not np.all(equal):
                distances = np.asarray(expected.get("unit_boundary_distances", []))
                eligible = (
                    np.any(distances <= atol, axis=1)
                    if distances.ndim == 2 and len(distances) == len(left)
                    else np.zeros(len(left), dtype=bool)
                )
                equal |= eligible
            if np.all(equal):
                continue
        elif key.endswith("::vertex_count"):
            if np.array_equal(left, right):
                continue
        else:
            left_finite, right_finite = np.isfinite(left), np.isfinite(right)
            if (
                not np.array_equal(left_finite, right_finite)
                or not np.array_equal(np.isnan(left), np.isnan(right))
                or not np.array_equal(np.isposinf(left), np.isposinf(right))
                or not np.array_equal(np.isneginf(left), np.isneginf(right))
            ):
                raise AssertionError(
                    f"{name} {key} nonfinite mask differs; max difference=inf at point 0"
                )
            mask = left_finite & right_finite
            delta = np.zeros(left.shape, dtype=float)
            delta[mask] = np.abs(left[mask] - right[mask])
            if key.endswith("::bounds") or key.endswith("::area"):
                zero_mismatch = (left == 0) != (right == 0)
                close = (
                    np.isclose(left, right, rtol=1e-6, atol=0.0, equal_nan=True) & ~zero_mismatch
                )
                close |= (left == 0) & (right == 0)
            else:
                close = (delta <= atol) | ~mask
            if np.all(close):
                continue
        try:
            delta = np.abs(left.astype(float) - right.astype(float))
            flat = int(np.nanargmax(delta)) if delta.size else 0
            maxdiff = float(np.nanmax(delta)) if delta.size else 0.0
            point = np.unravel_index(flat, delta.shape)[0] if delta.ndim else 0
        except (TypeError, ValueError):
            maxdiff, point = float("inf"), 0
        raise AssertionError(f"{name} {key} mismatch: max difference={maxdiff:g} at point {point}")
