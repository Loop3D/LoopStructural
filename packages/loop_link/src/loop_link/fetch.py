"""Decode viewer objects into scene objects.

The viewer sends an object as a JSON header and a set of named arrays. The
function :func:`decode_object` makes a scene object from them. The coordinates
of the result are world coordinates in float64.
"""

from __future__ import annotations

import itertools
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from loop_link.errors import InvalidInputError, UnsupportedObjectError
from loop_link.objects import (
    DEFAULT_CATEGORY_COLOURS,
    BlockModel,
    Categorical,
    Drillhole,
    Lines,
    Points,
    SceneObject,
    Surface,
    TetraMesh,
    Vectors,
)

__all__ = ["FetchResult", "decode_object"]

# The members of the property groups that make a Points object into a Vectors object.
_GROUP_MEMBERS = {
    "orientation": ("dip_direction", "dip"),
    "lineation": ("trend", "plunge"),
}

_SUPPORTED_TYPES = ("Points", "Curve", "Surface", "TetraMesh", "BlockModel", "Drillhole")


@dataclass(frozen=True, eq=False)
class FetchResult:
    """An object that was fetched from the viewer.

    Parameters
    ----------
    id : str
        The identifier of the object.
    object : SceneObject
        The decoded object, in world coordinates.
    object_type : str
        The type of the object in the viewer, for example ``"Curve"``.
    revision : int
        The revision of the object in the viewer.
    name : str or None
        The display name.
    parent : str or None
        The identifier of the parent group.
    origin : tuple of float or None
        The scene origin that was used to decode the object.
    crs : str or None
        The coordinate reference system of the scene.
    """

    id: str
    object: SceneObject
    object_type: str
    revision: int
    name: str | None
    parent: str | None
    origin: tuple[float, float, float] | None
    crs: str | None


def _warn_dropped(name: str, object_id: str, reason: str) -> None:
    """Tell the user that a property is not in the result."""
    warnings.warn(
        f"loop_link dropped the property '{name}' of object '{object_id}': {reason}",
        UserWarning,
        stacklevel=4,
    )


def _decode_text(desc: Mapping[str, Any], arrays: Mapping[str, np.ndarray]) -> np.ndarray:
    """Return a numpy unicode array from the bytes and the offsets of a text property."""
    data = np.asarray(arrays[desc["buffer"]], dtype=np.uint8).tobytes()
    offsets = np.asarray(arrays[desc["offsets"]]).astype(np.int64)
    items = [
        data[int(start) : int(end)].decode("utf-8") for start, end in itertools.pairwise(offsets)
    ]
    return np.array(items, dtype=str)


def _decode_referenced(
    desc: Mapping[str, Any], arrays: Mapping[str, np.ndarray], object_id: str
) -> Categorical:
    """Return a :class:`Categorical` from the codes and the value map of a property."""
    value_map = sorted(desc.get("value_map", []), key=lambda item: item["id"])
    ids = np.array([item["id"] for item in value_map], dtype=np.int64)
    names = [str(item["name"]) for item in value_map]
    colours = []
    for index, item in enumerate(value_map):
        colour = item.get("colour")
        if colour is None:
            colours.append(DEFAULT_CATEGORY_COLOURS[index % len(DEFAULT_CATEGORY_COLOURS)])
            continue
        values = tuple(int(c) for c in colour)
        colours.append(values if len(values) == 4 else (*values[:3], 255))
    raw = np.asarray(arrays[desc["buffer"]]).astype(np.int64)
    positions = np.searchsorted(ids, raw)
    safe = np.minimum(positions, max(len(ids) - 1, 0))
    if raw.size and (len(ids) == 0 or not np.all(ids[safe] == raw)):
        raise InvalidInputError(
            f"The property '{desc['name']}' of object '{object_id}' has a code that is not "
            "in its value map. Check the object in the viewer."
        )
    return Categorical(positions, names, colours)


def _decode_values(
    desc: Mapping[str, Any], arrays: Mapping[str, np.ndarray], object_id: str
) -> np.ndarray | Categorical:
    """Return the values of one data entry as a numpy array or a Categorical."""
    data_type = desc.get("data_type")
    if data_type == "referenced":
        return _decode_referenced(desc, arrays, object_id)
    if data_type == "text":
        return _decode_text(desc, arrays)
    values = np.asarray(arrays[desc["buffer"]])
    if data_type == "boolean":
        return values != 0
    if data_type == "integer":
        return values.astype(np.int32)
    if values.dtype == np.float32:
        return values
    return values.astype(np.float64)


def _route(
    desc: Mapping[str, Any],
    object_type: str,
    object_id: str,
    interval_set: str | None,
) -> str | None:
    """Return the target dictionary of a data entry, or None with a reason in a warning.

    The result is ``"properties"``, ``"cell_properties"`` or ``"interval_properties"``.
    """
    name = str(desc.get("name", desc.get("id")))
    association = desc.get("association")
    if association == "vertex":
        return "properties"
    if association == "cell":
        if object_type == "Drillhole":
            return "interval_properties"
        if object_type == "Points":
            raise UnsupportedObjectError(object_id, object_type)
        return "cell_properties"
    if association == "interval" and object_type == "Drillhole":
        log_set = desc.get("log_set")
        if log_set is None or log_set == interval_set:
            return "interval_properties"
        _warn_dropped(name, object_id, f"it belongs to the log set '{log_set}', not the first one")
        return None
    _warn_dropped(name, object_id, f"the association '{association}' is not supported")
    return None


def _group_pair(
    group: Mapping[str, Any],
    names: tuple[str, str],
    entries: Mapping[str, Mapping[str, Any]],
    arrays: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Return the two float64 arrays of a property group."""
    by_role = {member["role"]: member["data"] for member in group.get("members", [])}
    result = []
    for role in names:
        if role not in by_role or by_role[role] not in entries:
            raise InvalidInputError(
                f"The property group '{group.get('id')}' has no member '{role}'. "
                "Check the object in the viewer."
            )
        buffer = entries[by_role[role]]["buffer"]
        result.append(np.asarray(arrays[buffer]).astype(np.float64))
    return result[0], result[1]


def decode_object(
    header: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    origin: Sequence[float] | None,
) -> SceneObject:
    """Decode the header and the arrays of a viewer object into a scene object.

    Parameters
    ----------
    header : mapping
        The ``scene.upsert`` header of the object.
    arrays : mapping of str to numpy.ndarray
        The buffers of the object, by name.
    origin : sequence of float or None
        The scene origin. The function adds it to the coordinates. None means
        (0, 0, 0).

    Returns
    -------
    SceneObject
        The object, with float64 world coordinates.

    Raises
    ------
    UnsupportedObjectError
        If the object type is not supported, or a cell property is on a Points object.
    InvalidInputError
        If a code of a referenced property is not in its value map.

    Warns
    -----
    UserWarning
        One time for each property that the function drops.
    """
    object_id = str(header.get("id"))
    object_type = str(header.get("object_type"))
    if object_type not in _SUPPORTED_TYPES:
        raise UnsupportedObjectError(object_id, object_type)
    shift = np.zeros(3) if origin is None else np.asarray(origin, dtype=np.float64)
    geometry = header.get("geometry", {})
    params = header.get("params", {})

    def geo(role: str) -> np.ndarray:
        return np.asarray(arrays[geometry[role]])

    entries = {str(desc["id"]): desc for desc in header.get("data", [])}
    groups = [
        group
        for group in header.get("property_groups", [])
        if object_type == "Points" and group.get("role") in _GROUP_MEMBERS
    ]
    member_ids = {member["data"] for group in groups for member in group.get("members", [])}

    interval_set = None
    for log_set in header.get("log_sets", []):
        if log_set.get("kind") == "interval":
            interval_set = log_set.get("name")
            break

    targets: dict[str, dict[str, Any]] = {
        "properties": {},
        "cell_properties": {},
        "interval_properties": {},
    }
    for data_id, desc in entries.items():
        if data_id in member_ids:
            continue
        target = _route(desc, object_type, object_id, interval_set)
        if target is not None:
            targets[target][str(desc["name"])] = _decode_values(desc, arrays, object_id)

    if object_type == "Drillhole":
        end_of_hole = params.get("end_of_hole")
        return Drillhole(
            geo("collar").astype(np.float64) + shift,
            geo("survey").astype(np.float64),
            geo("intervals").astype(np.float64),
            end_of_hole=None if end_of_hole is None else float(end_of_hole[0]),
            properties=targets["properties"],
            interval_properties=targets["interval_properties"],
        )
    if object_type == "BlockModel":
        return BlockModel(
            np.asarray(params["origin"], dtype=np.float64) + shift,
            params["cell_sizes"],
            params["cells"],
            properties=targets["properties"],
            cell_properties=targets["cell_properties"],
        )

    vertices = geo("vertices").astype(np.float64) + shift
    if object_type == "Curve":
        return Lines(vertices, geo("segments"), targets["properties"], targets["cell_properties"])
    if object_type == "Surface":
        return Surface(
            vertices, geo("triangles"), targets["properties"], targets["cell_properties"]
        )
    if object_type == "TetraMesh":
        return TetraMesh(
            vertices, geo("tetrahedra"), targets["properties"], targets["cell_properties"]
        )
    if not groups:
        return Points(vertices, targets["properties"])
    pairs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for group in groups:
        role = group["role"]
        pairs[role] = _group_pair(group, _GROUP_MEMBERS[role], entries, arrays)
    return Vectors(
        vertices,
        orientation=pairs.get("orientation"),
        lineation=pairs.get("lineation"),
        properties=targets["properties"],
    )
