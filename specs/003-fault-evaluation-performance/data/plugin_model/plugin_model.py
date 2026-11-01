"""Build the plugin model of spec 003 (research.md) from the files in this folder.

The QGIS plugin built this model. The files are:

- ``recipe.json``: the box, the stratigraphic column, and the LoopStructural
  calls that the plugin made, in order (recorded from the plugin on
  2026-11-01).
- ``<feature>.csv``: the data of each feature.
- ``fields.npz`` and ``fields.json``: the solved field of each interpolator
  (one for each coordinate of a fault frame), and its grid.
- ``qgis_check.npz``: values of the model in QGIS at 5000 fixed points, to
  check a build.

``build_model()`` gives the model with the saved fields, so nothing is solved
and the result does not depend on the solver. ``build_model(solve=True)``
solves the model again from the data.

Run this file to check that a build gives the values of the model in QGIS.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.geometry import BoundingBox
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn

HERE = Path(__file__).parent


def _builders(feature: object) -> list:
    builder = feature.builder
    return list(getattr(builder, "builders", None) or [builder])


def _load_fields(feature: object, fields: np.lib.npyio.NpzFile, info: dict) -> None:
    """Put the saved field into each interpolator of ``feature``, and mark it as built."""
    builders = _builders(feature)
    for index, builder in enumerate(builders):
        key = f"{feature.name}__{index}" if len(builders) > 1 else feature.name
        interpolator = builder.interpolator
        # The build makes the grid of a fault larger (the fault buffer), so set
        # the saved grid. The setters change each other, so set the values.
        support = interpolator.support
        support._origin = np.array(info[key]["origin"], dtype=float)
        support._step_vector = np.array(info[key]["step_vector"], dtype=float)
        support._nsteps = np.array(info[key]["nsteps"], dtype=np.int64)
        support.onGeometryChange()
        if int(interpolator.support.n_nodes) != info[key]["n_nodes"]:
            raise ValueError(
                f"The grid of {key} has {interpolator.support.n_nodes} nodes, "
                f"the saved field has {info[key]['n_nodes']}."
            )
        interpolator.c = np.array(fields[key], dtype=float)
        interpolator.up_to_date = True
        builder._up_to_date = True
    feature.builder._up_to_date = True


def _copy_faults(feature: object) -> None:
    """Give the feature the faults of its builder, as ``BaseBuilder.update`` does.

    Do this after all the calls: ``add_unconformity`` must see the feature
    before its first build, as in the plugin.
    """
    for builder in _builders(feature):
        builder._feature.faults = builder.faults


def build_model(solve: bool = False, solver: str | None = None) -> GeologicalModel:
    """Build the plugin model.

    Parameters
    ----------
    solve : bool
        False: use the saved fields (no solve). True: solve again.
    solver : str, optional
        The solver for ``solve=True`` (for example ``"cg"`` or ``"native"``).
    """
    recipe = json.loads((HERE / "recipe.json").read_text())
    info = json.loads((HERE / "fields.json").read_text())
    fields = np.load(HERE / "fields.npz")
    # The plugin gives the model a BoundingBox, so its local origin stays at 0
    # (world coordinates). GeologicalModel(origin, maximum) would move it.
    model = GeologicalModel(BoundingBox.from_dict(recipe["bounding_box"]))
    model.stratigraphic_column = StratigraphicColumn.from_dict(recipe["stratigraphic_column"])
    loaded = []
    for step in recipe["steps"]:
        call = step["call"]
        if call in ("create_and_add_fault", "create_and_add_foliation"):
            kwargs = dict(step["kwargs"])
            if "faults" in kwargs:
                kwargs["faults"] = [model.get_feature_by_name(n) for n in kwargs["faults"]]
            # the grid of the model in QGIS (the user can change it in the plugin)
            key = f"{step['name']}__0" if call == "create_and_add_fault" else step["name"]
            kwargs["nelements"] = info[key]["nelements"]
            if solve and solver is not None:
                kwargs["solver"] = solver
            data = pd.read_csv(HERE / step["data"], dtype={"feature_name": str})
            feature = getattr(model, call)(step["name"], data=data, **kwargs)
            if not solve:
                _load_fields(feature, fields, info)
                loaded.append(feature)
        elif call == "add_abutting_fault":
            fault = model.get_feature_by_name(step["fault"])
            fault.add_abutting_fault(
                model.get_feature_by_name(step["abutting_fault"]), positive=step["positive"]
            )
        elif call in ("add_unconformity", "add_onlap_unconformity"):
            getattr(model, call)(model.get_feature_by_name(step["feature"]), step["value"])
        else:
            raise ValueError(f"Unknown call {call}")
    for feature in loaded:
        _copy_faults(feature)
    if solve:
        model.update()
    return model


def check(model: GeologicalModel | None = None) -> dict:
    """Return the largest difference from the values of the model in QGIS."""
    model = model or build_model()
    reference = np.load(HERE / "qgis_check.npz")
    points = reference["points"]
    result = {}
    group = model.get_feature_by_name("Group_0")
    result["Group_0"] = group.evaluate_value(points)
    for name in ("7439", "2997", "3496", "3498", "12647", "12658"):
        result[f"{name}__0"] = model.get_feature_by_name(name)[0].evaluate_value(points)
    result["restored_Group_0"] = group._apply_faults(points.copy())
    differences = {}
    for key, value in result.items():
        expected = reference[key]
        same_nan = np.array_equal(np.isnan(expected), np.isnan(value))
        both = ~np.isnan(expected) & ~np.isnan(value)
        largest = float(np.max(np.abs(expected[both] - value[both]))) if both.any() else 0.0
        differences[key] = (same_nan, largest)
    return differences


if __name__ == "__main__":
    for key, (same_nan, largest) in check().items():
        print(f"{key}: same NaN: {same_nan}, largest difference: {largest:.3g}")
