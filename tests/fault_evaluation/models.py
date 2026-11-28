"""Deterministic fixed-field models used by the fault evaluation suite."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.geometry import BoundingBox
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn
from LoopStructural.modelling.features.fold.fold_function import (
    FoldRotationType,
    get_fold_rotation_profile,
)

_HERE = Path(__file__).parent
_DATA = _HERE / "data"
_RECIPES = json.loads((_DATA / "recipes.json").read_text(encoding="utf8"))
_SYNTHETIC = _RECIPES["synthetic"]
_SPECIAL = ("stopped_fault", "cutting_faults", "onlap_domain", "faulted_fold")
_FIELD_FILE = _DATA / "fields.npz"


def model_names() -> list[str]:
    """Return the fixed catalogue in baseline order."""
    return [
        f"f{faults}_u{unconformities}"
        for faults in _SYNTHETIC["faults"]
        for unconformities in _SYNTHETIC["unconformities"]
    ] + [*_SPECIAL, "plugin"]


def _rows(fault_count: int, unconformity_count: int) -> pd.DataFrame:
    rng = np.random.default_rng(_SYNTHETIC["seed"])
    rows: list[dict[str, Any]] = []
    for i in range(fault_count):
        x = 150 + 700 * (i + 0.5) / max(fault_count, 1)
        name = f"fault_{i}"
        for y in (100.0, 500.0, 900.0):
            rows.append({"X": x, "Y": y, "Z": -500.0, "val": 0.0, "feature_name": name})
        rows.append(
            {
                "X": x,
                "Y": 500.0,
                "Z": -500.0,
                "gx": 1.0,
                "gy": 0.2 * i,
                "gz": 0.0,
                "feature_name": name,
            }
        )
    for k in range(unconformity_count + 1):
        name, tilt = f"group_{k}", 0.05 * k
        for value in (0.0, 50.0, 100.0):
            xy = rng.random((30, 2)) * 1000
            z = -900 + 150 * k + value + tilt * xy[:, 0]
            for (x, y), zz in zip(xy, z):
                rows.append({"X": x, "Y": y, "Z": zz, "val": value, "feature_name": name})
        rows.append(
            {
                "X": 500.0,
                "Y": 500.0,
                "Z": -900 + 150 * k,
                "gx": -tilt,
                "gy": 0.0,
                "gz": 1.0,
                "feature_name": name,
            }
        )
    return pd.DataFrame(rows)


def _builders(feature: Any) -> list[Any]:
    """Return the real builders owned by a feature or its coordinates once."""
    result: list[Any] = []
    seen: set[int] = set()
    for coordinate in list(getattr(feature, "features", None) or []) + [feature]:
        owner = getattr(coordinate, "builder", None)
        for builder in list(getattr(owner, "builders", None) or [owner]):
            if builder is not None and id(builder) not in seen:
                seen.add(id(builder))
                result.append(builder)
    return result


def _load_plugin() -> GeologicalModel:
    path = (
        _HERE.parent.parent
        / "specs"
        / "003-fault-evaluation-performance"
        / "data"
        / "plugin_model"
        / "plugin_model.py"
    )
    spec = importlib.util.spec_from_file_location("_fault_evaluation_plugin_model", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load supplied plugin builder at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_model()


def _build_recipe(name: str, *, solve: bool = False) -> GeologicalModel:
    if name == "plugin":
        return _load_plugin()
    special = name in _SPECIAL
    if special:
        fault_count = int(_RECIPES["special"][name].get("faults", 0))
        unconformity_count = int(_RECIPES["special"][name].get("unconformities", 0))
    else:
        fault_count, unconformity_count = int(name[1]), int(name.rsplit("_u", 1)[1])
    origin, maximum = np.array(_SYNTHETIC["origin"]), np.array(_SYNTHETIC["maximum"])
    if name == "faulted_fold":
        _, bounds = load_noddy_single_fold()
        origin, maximum = bounds[0].copy(), bounds[1].copy()
    if name in ("cutting_faults", "faulted_fold"):
        # Keep model coordinates in world space for nested fault-frame builds
        # and the folded Noddy frame. The builder validates the frame by
        # restoring a buffered local grid; a translated frame would translate
        # that grid a second time.
        box = BoundingBox(origin, maximum)
        box.set_local_transform(
            local_origin=np.asarray(_RECIPES["special"][name]["local_origin"], dtype=float)
        )
        model = GeologicalModel(box)
    else:
        model = GeologicalModel(origin, maximum)
    data = _rows(fault_count, unconformity_count)
    if name == "faulted_fold":
        fold_data, _ = load_noddy_single_fold()
        # Preserve ``coord``: the shipped dataset uses it for two s1 frame
        # constraints, while s0 values are supplied through other columns.
        fold_data["feature_name"] = fold_data["feature_name"].replace(
            {"s1": "fold_frame", "s0": "folded_unit"}
        )
        # The supplied dataset has validated coordinate and orientation constraints.
        # Add a fault trace and dip measurement inside its recorded box.
        fault_x = float(_RECIPES["special"][name]["fault_x"])
        fault_rows = [
            {"X": fault_x, "Y": y, "Z": 7500.0, "val": 0.0, "coord": 0, "feature_name": "fault_0"}
            for y in (1000.0, 3500.0, 6000.0)
        ]
        fault_rows.append(
            {
                "X": fault_x,
                "Y": 3500.0,
                "Z": 7500.0,
                "gx": 1.0,
                "gy": 0.0,
                "gz": 0.0,
                "coord": 0,
                "feature_name": "fault_0",
            }
        )
        group_rows = []
        rng = np.random.default_rng(31)
        for value in (0.0, 50.0, 100.0):
            xy = rng.uniform(origin[:2], maximum[:2], (40, 2))
            z = origin[2] + 500 + value + 0.05 * (xy[:, 0] - origin[0])
            group_rows.extend(
                {"X": x, "Y": y, "Z": zz, "val": value, "feature_name": "group_0"}
                for (x, y), zz in zip(xy, z)
            )
        data = pd.concat(
            (fold_data, pd.DataFrame(fault_rows + group_rows)), ignore_index=True, sort=False
        )
    if name == "onlap_domain":
        rng = np.random.default_rng(29)
        rows = []
        for value in (0.0, 50.0, 100.0):
            xy = rng.uniform(0, 1000, (40, 2))
            rows.extend(
                {
                    "X": x,
                    "Y": y,
                    "Z": -800 + value + 0.12 * x,
                    "val": value,
                    "feature_name": "domain_surface",
                }
                for x, y in xy
            )
        data = pd.concat((data, pd.DataFrame(rows)), ignore_index=True)
    model.data = data
    fault_names = [f"fault_{i}" for i in range(fault_count)]
    groups = []
    if name == "onlap_domain":
        # The domain-fault factory reads model.data directly. Prepare the raw
        # measurements first so absent normal and vector columns are NaN.
        # Add it before fault_0: its mask then has no dependency on a fault
        # that the domain surface would later crop.
        model.data = model.prepare_data(model.data)
        model.create_and_add_domain_fault("domain_surface", nelements=_SYNTHETIC["nelements"])
    if name == "stopped_fault":
        younger = model.create_and_add_foliation("group_0", nelements=_SYNTHETIC["nelements"])
        groups.append(younger)
        horizon = model.add_unconformity(younger, 100.0)
        if fault_names:
            model.create_and_add_fault(
                fault_names[0], _SYNTHETIC["fault_displacement"], nelements=_SYNTHETIC["nelements"]
            )
        groups.append(model.create_and_add_foliation("group_1", nelements=_SYNTHETIC["nelements"]))
    else:
        for i, fault in enumerate(fault_names):
            kwargs = {"nelements": _SYNTHETIC["nelements"]}
            if name == "cutting_faults" and i:
                kwargs["faults"] = [model.get_feature_by_name("fault_0")]
            if name == "faulted_fold":
                kwargs["fault_center"] = _RECIPES["special"][name]["fault_center"]
            model.create_and_add_fault(fault, _SYNTHETIC["fault_displacement"], **kwargs)
        for k in range(unconformity_count + 1):
            feature = model.create_and_add_foliation(
                f"group_{k}", nelements=_SYNTHETIC["nelements"]
            )
            groups.append(feature)
            if k < unconformity_count:
                if name == "onlap_domain":
                    model.add_onlap_unconformity(feature, 100.0)
                else:
                    model.add_unconformity(feature, 100.0)
    if name == "faulted_fold":
        frame = model.create_and_add_fold_frame("fold_frame", nelements=_SYNTHETIC["nelements"])
        folded = model.create_and_add_folded_foliation(
            "folded_unit",
            fold_frame=frame,
            nelements=_SYNTHETIC["nelements"],
            fold_axis=_RECIPES["special"][name]["fold_axis"],
            limb_wl=_RECIPES["special"][name]["limb_wavelength"],
        )
    # Declare two real scalar intervals for each feature group. Add them in
    # reverse group order because StratigraphicColumn.get_groups() reads the
    # column from the top back toward the oldest group.
    column = StratigraphicColumn()
    column.clear(basement=False)
    stratigraphy = _RECIPES["synthetic"]["stratigraphic_column"]
    unit_intervals: list[tuple[Any, float, float]] = []
    for k in reversed(range(unconformity_count + 1)):
        if k == 0 and unconformity_count:
            bounds_key = "onlap_youngest" if name == "onlap_domain" else "erosional_youngest"
            bounds = stratigraphy[bounds_key]
        elif k == 0:
            bounds = stratigraphy["default"]
        else:
            low, high = stratigraphy["older_outer"]
            contact = stratigraphy["older_contact_origin"] + (
                stratigraphy["older_contact_per_index"] * k
            )
            bounds = [low, contact, high]
        low, contact, high = map(float, bounds)
        # Reverse insertion order so each group's `units` list is low then high.
        high_unit = column.add_unit(f"unit_{k}_high", thickness=high - contact)
        low_unit = column.add_unit(f"unit_{k}_low", thickness=contact - low)
        unit_intervals.extend(((high_unit, contact, high), (low_unit, low, contact)))
        if k > 0:
            column.add_unconformity(f"unconformity_{k}")
    # add_unit/add_unconformity recalculate cumulative thickness, so install
    # the declared scalar bounds only after the order is complete.
    for unit, minimum, maximum_value in unit_intervals:
        unit.min_value = minimum
        unit.max_value = maximum_value
    column.group_mapping.update({f"Group_{k}": f"group_{k}" for k in range(unconformity_count + 1)})
    model.stratigraphic_column = column
    # Match the plugin loading contract: attach each builder's fault list only
    # after all relationships have been constructed. The normal path never solves.
    for feature in model.features:
        for builder in _builders(feature):
            if hasattr(builder, "_feature") and hasattr(builder, "faults"):
                builder._feature.faults = builder.faults
    if solve:
        if name == "faulted_fold":
            recipe = _RECIPES["special"][name]
            frame.builder.up_to_date()
            folded.builder.set_fold_axis()
            angles, coordinates = folded.builder.calculate_fold_limb_rotation_angle()
            profile = get_fold_rotation_profile(
                folded.builder.limb_profile_type, angles, coordinates
            )
            profile.fit(
                params={
                    "wavelength": recipe["limb_wavelength"],
                    "calculate_wavelength": recipe["calculate_wavelength"],
                }
            )
            folded.fold.fold_limb_rotation = profile
            profile.add_observer(folded.builder)
        model.update()
    return model


def _interpolated_features(model: GeologicalModel) -> list[tuple[Any, Any]]:
    """Return each real feature-coordinate field and its owning builder once."""
    result = []
    seen: set[int] = set()
    for feature in model.features:
        candidates = list(getattr(feature, "features", None) or []) + [feature]
        for coordinate in candidates:
            builder = getattr(coordinate, "builder", None)
            if builder is None or not hasattr(builder, "interpolator") or id(builder) in seen:
                continue
            if not hasattr(builder.interpolator, "support"):
                continue
            seen.add(id(builder))
            result.append((coordinate, builder))
    return result


def _field_key(name: str, feature: Any) -> str:
    return f"{name}::{feature.name}"


def _apply_fields(model: GeologicalModel, name: str, archive: Any) -> None:
    metadata = json.loads(str(archive["__metadata__"]))
    for feature, builder in _interpolated_features(model):
        key = _field_key(name, feature)
        if key not in metadata:
            raise KeyError(f"Stored field missing for {key}")
        info = metadata[key]
        support = builder.interpolator.support
        if info["support_type"] != type(support).__name__:
            raise TypeError(f"Stored support type differs for {key}")
        support._origin = np.asarray(info["origin"], dtype=float)
        support._step_vector = np.asarray(info["step_vector"], dtype=float)
        support._nsteps = np.asarray(info["nsteps"], dtype=np.int64)
        if hasattr(support, "rotation_xy"):
            support.rotation_xy = np.asarray(info["rotation_xy"], dtype=float)
        support.onGeometryChange()
        builder.interpolator.c = archive[key]
        builder.interpolator.up_to_date = info["interpolator_up_to_date"]
        builder._up_to_date = info["builder_up_to_date"]
    for feature in model.features:
        if getattr(feature, "builder", None) is not None:
            feature.builder._up_to_date = True
    fold_key = f"{name}::__fold_events__"
    if fold_key in metadata:
        profiles = {
            profile.__name__: profile
            for profile in (
                FoldRotationType.FOURIER_SERIES.value,
                FoldRotationType.TRIGONOMETRIC.value,
            )
        }
        for event_info in metadata[fold_key]["events"]:
            feature = model.get_feature_by_name(event_info["feature"])
            fold = feature.fold
            fold.fold_axis = (
                None if event_info["fold_axis"] is None else np.asarray(event_info["fold_axis"])
            )
            fold.invert_norm = bool(event_info["invert_norm"])
            for attribute in ("fold_axis_rotation", "fold_limb_rotation"):
                profile_info = event_info.get(attribute)
                if profile_info is None:
                    setattr(fold, attribute, None)
                    continue
                profile_type = profiles[profile_info["type"]]
                rotation_type = next(
                    item for item in FoldRotationType if item.value is profile_type
                )
                profile = get_fold_rotation_profile(
                    rotation_type,
                    np.asarray(profile_info["rotation_angle"], dtype=float),
                    np.asarray(profile_info["fold_frame_coordinate"], dtype=float),
                )
                profile.update_params(list(profile_info["params"].values()))
                profile.evaluation_points = np.asarray(
                    profile_info["evaluation_points"], dtype=float
                )
                setattr(fold, attribute, profile)
                profile.add_observer(feature.builder)


def _save_solved_fields(path: Path = _FIELD_FILE) -> None:
    """Internal one-time helper: solve each recipe and freeze its grid and coefficients."""
    arrays: dict[str, Any] = {}
    metadata: dict[str, Any] = {}
    for name in model_names():
        model = _build_recipe(name, solve=True)
        if name == "plugin":
            # The supplied plugin loader consumes its original, immutable archive.
            continue
        for feature, builder in _interpolated_features(model):
            interp = builder.interpolator
            key = _field_key(name, feature)
            arrays[key] = np.asarray(interp.c)
            metadata[key] = {
                "support_type": type(interp.support).__name__,
                "origin": np.asarray(interp.support.origin).tolist(),
                "step_vector": np.asarray(interp.support.step_vector).tolist(),
                "nsteps": np.asarray(interp.support.nsteps).tolist(),
                "rotation_xy": np.asarray(
                    getattr(interp.support, "rotation_xy", np.eye(3))
                ).tolist(),
                "builder_up_to_date": bool(builder._up_to_date),
                "interpolator_up_to_date": bool(interp.up_to_date),
            }
        if name == "faulted_fold":
            folds = []
            for feature in model.features:
                fold = getattr(feature, "fold", None)
                if fold is not None:

                    def profile_state(profile: Any) -> dict[str, Any] | None:
                        if profile is None:
                            return None
                        return {
                            "type": type(profile).__name__,
                            "rotation_angle": np.asarray(profile.rotation_angle).tolist(),
                            "fold_frame_coordinate": np.asarray(
                                profile.fold_frame_coordinate
                            ).tolist(),
                            "params": {key: float(value) for key, value in profile.params.items()},
                            "evaluation_points": np.asarray(profile.evaluation_points).tolist(),
                        }

                    folds.append(
                        {
                            "feature": feature.name,
                            "frame": fold.foldframe.name,
                            "name": fold.name,
                            "fold_axis": (
                                None
                                if fold.fold_axis is None
                                else np.asarray(fold.fold_axis).tolist()
                            ),
                            "invert_norm": bool(fold.invert_norm),
                            "fold_axis_rotation": profile_state(fold.fold_axis_rotation),
                            "fold_limb_rotation": profile_state(fold.fold_limb_rotation),
                        }
                    )
            metadata[f"{name}::__fold_events__"] = {"events": folds}
    arrays["__metadata__"] = np.asarray(json.dumps(metadata))
    np.savez_compressed(path, **arrays)


def build_model(name: str) -> GeologicalModel:
    """Create a model and load its saved fields without running any solver."""
    if name not in model_names():
        raise ValueError(f"Unknown fault evaluation model {name!r}")
    model = _build_recipe(name)
    if name == "plugin":
        return model
    if not _FIELD_FILE.exists():
        raise FileNotFoundError(f"Golden field archive is missing: {_FIELD_FILE}")
    with np.load(_FIELD_FILE) as archive:
        _apply_fields(model, name, archive)
    return model


def sample_points(name: str, count: int, seed: int = 1) -> np.ndarray:
    """Return deterministic uniform points within the named model box."""
    model = build_model(name)
    rng = np.random.default_rng(seed)
    box = model.bounding_box
    return rng.random((count, 3)) * (box.maximum - box.origin) + box.origin
