"""Capture and compare deterministic fault-evaluation outputs."""

from __future__ import annotations

from typing import Any
from unittest.mock import patch

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


def _surface_sample_capture(feature: Any, points: np.ndarray):
    """Return a scoped evaluator wrapper for known pointwise geological fields."""
    try:
        from LoopStructural.modelling.features._base_geological_feature import BaseFeature
        from LoopStructural.modelling.features._geological_feature import GeologicalFeature
        from LoopStructural.modelling.features._unconformity_feature import UnconformityFeature
        from loop_common.supports._3d_structured_grid import StructuredGrid as SupportGrid
        from loop_interpolation._discrete_interpolator import DiscreteInterpolator
        from loop_interpolation._geological_interpolator import GeologicalInterpolator
    except ImportError:
        return None

    evaluator = getattr(feature, "evaluate_value", None)
    interpolator = getattr(feature, "interpolator", None)
    if (
        not isinstance(feature, GeologicalFeature)
        or getattr(evaluator, "__func__", None) is not GeologicalFeature.evaluate_value
        or getattr(getattr(feature, "_calculate_mask", None), "__func__", None)
        is not BaseFeature._calculate_mask
        or getattr(getattr(feature, "_apply_faults", None), "__func__", None)
        is not BaseFeature._apply_faults
        or not isinstance(interpolator, DiscreteInterpolator)
        or getattr(getattr(interpolator, "_evaluate_value_local", None), "__func__", None)
        is not DiscreteInterpolator._evaluate_value_local
        or getattr(getattr(interpolator, "evaluate_value", None), "__func__", None)
        is not GeologicalInterpolator.evaluate_value
        or type(getattr(interpolator, "support", None)) is not SupportGrid
        or any(
            type(region) is not UnconformityFeature
            or getattr(getattr(region, "evaluate_value", None), "__func__", None)
            is not GeologicalFeature.evaluate_value
            for region in (getattr(feature, "regions", None) or [])
        )
    ):
        return None

    original = evaluator
    recorded: dict[str, np.ndarray] = {}
    depth = 0

    def record(pos, *args, **kwargs):
        nonlocal depth
        outermost = depth == 0
        depth += 1
        try:
            values = original(pos, *args, **kwargs)
            ignore_regions = kwargs.get("ignore_regions", args[0] if args else False)
            fillnan = kwargs.get("fillnan", args[1] if len(args) > 1 else None)
            array = np.asarray(pos)
            if (
                outermost
                and not recorded
                and points.shape == (8000, 3)
                and array.shape == points.shape
                and np.array_equal(array, points)
                and ignore_regions is False
                and fillnan is None
            ):
                value_array = np.asarray(values)
                if value_array.shape == (8000,):
                    recorded["values"] = value_array.copy()
            return values
        finally:
            depth -= 1

    return patch.object(feature, "evaluate_value", record), recorded


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
        surface_points = box.regular_grid(local=False, order="C")
        surface_names = [f"golden-level:{level!r}" for level in group_levels]
        scoped_capture = _surface_sample_capture(feature, surface_points)
        if scoped_capture is None:
            meshes = list(feature.surfaces(group_levels, box, name=surface_names))
            recorded_sample = None
        else:
            evaluator_patch, recorded_sample = scoped_capture
            with evaluator_patch:
                meshes = list(feature.surfaces(group_levels, box, name=surface_names))

        names_to_levels = dict(zip(surface_names, group_levels))
        meshes_by_level: dict[float, Any] = {}
        all_metadata_absent = bool(meshes) and all(
            not hasattr(mesh, "values") and not hasattr(mesh, "name") for mesh in meshes
        )
        if all_metadata_absent:
            if len(meshes) != len(group_levels):
                raise ValueError(
                    f"Ambiguous surface count for {group.name}: {len(meshes)} != {len(group_levels)}"
                )
            meshes_by_level = dict(zip(group_levels, meshes))
        else:
            for mesh in meshes:
                vertices = np.asarray(getattr(mesh, "vertices", np.empty((0, 3))))
                mesh_name = getattr(mesh, "name", None)
                named_level = names_to_levels.get(mesh_name)
                has_values = hasattr(mesh, "values")
                values = np.asarray(getattr(mesh, "values", []), dtype=float).reshape(-1)

                if values.size:
                    if not np.all(np.isfinite(values)) or not np.all(values == values[0]):
                        raise ValueError(
                            f"Non-constant or non-finite surface value for {group.name}"
                        )
                    mesh_level = float(values[0])
                    if mesh_level not in group_levels:
                        raise ValueError(f"Unexpected surface value for {group.name}: {mesh_level}")
                    if mesh_name is not None and (named_level is None or named_level != mesh_level):
                        raise ValueError(f"Surface name/value mismatch for {group.name}")
                elif len(vertices) == 0 and has_values and named_level is not None:
                    # Actual Surface removes NaN vertices during construction,
                    # leaving an empty values array. The explicit requested
                    # name is the only reliable level identity in that case.
                    mesh_level = named_level
                else:
                    raise ValueError(f"Unrecognized surface level for {group.name}")

                if mesh_level in meshes_by_level:
                    raise ValueError(f"Duplicate surface value for {group.name}: {mesh_level}")
                meshes_by_level[mesh_level] = mesh
        grid_points = (
            surface_points if recorded_sample is not None and "values" in recorded_sample else None
        )
        finite_min = finite_max = None
        visible = None
        mask_checked = False
        if recorded_sample is not None and "values" in recorded_sample:
            grid_points = surface_points
            original_mask = feature._calculate_mask(grid_points)
            diagnostic_values = recorded_sample["values"]
            if diagnostic_values.shape == (len(grid_points),):
                diagnostic_finite = diagnostic_values[
                    original_mask & np.isfinite(diagnostic_values)
                ]
                if diagnostic_finite.size:
                    finite_min = float(diagnostic_finite.min())
                    finite_max = float(diagnostic_finite.max())
                visible = original_mask
                mask_checked = True
        for level in group_levels:
            mesh = meshes_by_level.get(level)
            if mesh is None:
                vertices = np.empty((0, 3))
                triangles = np.empty((0, 3), dtype=int)
            else:
                vertices = np.asarray(getattr(mesh, "vertices", np.empty((0, 3))))
                triangles = np.asarray(getattr(mesh, "triangles", np.empty((0, 3), dtype=int)))
            if len(vertices):
                _store_surface(result, feature.name, level, vertices, triangles)
                continue

            # Reuse one field sample and its mask diagnostics for all empty levels.
            if finite_min is None and not mask_checked:
                if grid_points is None:
                    grid_points = box.regular_grid(local=False, order="F")
                samples = np.asarray(feature.evaluate_value(grid_points))
                finite = samples[np.isfinite(samples)]
                if finite.size:
                    finite_min, finite_max = float(finite.min()), float(finite.max())
            outside_range = finite_min is not None and (level < finite_min or level > finite_max)
            fully_eroded = False
            if not outside_range and not mask_checked:
                mask_fn = getattr(feature, "_calculate_mask", None)
                if grid_points is None:
                    grid_points = box.regular_grid(local=False, order="F")
                visible = mask_fn(grid_points) if callable(mask_fn) else None
                mask_checked = True
            if not outside_range and visible is not None:
                fully_eroded = not np.any(visible)
            if not fully_eroded and not outside_range:
                raise ValueError(f"Empty surface for {group.name} at {level:g}")
            _store_surface(result, feature.name, level, None, None)

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
