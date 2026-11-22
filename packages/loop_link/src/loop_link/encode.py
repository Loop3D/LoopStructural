"""Encode scene objects into protocol frames.

The public functions return ``bytes``. Each frame has a JSON header and a list of
named binary buffers. See :mod:`loop_link.protocol` for the frame layout.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from loop_link import protocol
from loop_link.errors import InvalidInputError, ObjectTooLargeError
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

__all__ = [
    "cancel_frame",
    "clear_frame",
    "command_frame",
    "commands_list_frame",
    "encode_object",
    "fetch_frame",
    "group_frame",
    "heartbeat_frame",
    "hello_frame",
    "remove_frame",
    "selection_get_frame",
    "selection_set_frame",
]

_I32_MIN = -(2**31)
_I32_MAX = 2**31 - 1

# The name of the log set that holds the interval properties of a Drillhole.
_INTERVAL_LOG_SET = "intervals"

Buffers = list[tuple[str, np.ndarray]]


def hello_frame(crs: str, origin: Sequence[float], version: str) -> bytes:
    """Return the ``hello`` frame.

    Parameters
    ----------
    crs : str
        The coordinate reference system, for example ``"EPSG:28350"``.
    origin : sequence of float
        The scene origin, three values.
    version : str
        The version of loop_link.

    Returns
    -------
    bytes
        The frame.
    """
    header = {
        "type": "hello",
        "application": {"name": "loop_link", "version": version},
        "crs": crs,
        "origin": [float(v) for v in origin],
    }
    return protocol.encode_frame(header)


def heartbeat_frame(seq: int) -> bytes:
    """Return a ``heartbeat`` frame.

    Parameters
    ----------
    seq : int
        The sequence number.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame({"type": "heartbeat", "seq": seq})


def remove_frame(id: str) -> bytes:
    """Return a ``scene.remove`` frame.

    Parameters
    ----------
    id : str
        The identifier of the object to remove.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame({"type": "scene.remove", "id": id})


def clear_frame() -> bytes:
    """Return a ``scene.clear`` frame.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame({"type": "scene.clear"})


def selection_set_frame(
    request: int,
    object_id: str,
    cells: np.ndarray | None,
    vertices: np.ndarray | None,
    revision: int | None,
) -> bytes:
    """Return a ``selection.set`` frame.

    An empty or None array adds no entry and no buffer. The function does not
    check the index range.

    Parameters
    ----------
    request : int
        The request number.
    object_id : str
        The identifier of the object.
    cells : numpy.ndarray or None
        The cell indices.
    vertices : numpy.ndarray or None
        The vertex (node) indices.
    revision : int or None
        The revision of the object. If it is None, the header has no ``revision`` key.

    Returns
    -------
    bytes
        The frame.
    """
    header: dict[str, Any] = {"type": "selection.set", "request": request, "object": object_id}
    if revision is not None:
        header["revision"] = revision
    entries: list[dict[str, str]] = []
    buffers: Buffers = []
    for kind, name, values in (("cells", "cells", cells), ("nodes", "nodes", vertices)):
        if values is None:
            continue
        array = np.ascontiguousarray(np.asarray(values, dtype=np.uint32))
        if array.size == 0:
            continue
        entries.append({"kind": kind, "indices": name})
        buffers.append((name, array))
    header["entries"] = entries
    return protocol.encode_frame(header, buffers, version=protocol.TOOL_PROTOCOL_VERSION)


def selection_get_frame(request: int) -> bytes:
    """Return a ``selection.get`` frame.

    Parameters
    ----------
    request : int
        The request number.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame(
        {"type": "selection.get", "request": request},
        version=protocol.TOOL_PROTOCOL_VERSION,
    )


def fetch_frame(request: int, id: str) -> bytes:
    """Return an ``object.fetch`` frame.

    Parameters
    ----------
    request : int
        The request number.
    id : str
        The identifier of the object.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame(
        {"type": "object.fetch", "request": request, "id": id},
        version=protocol.TOOL_PROTOCOL_VERSION,
    )


def cancel_frame(request: int) -> bytes:
    """Return a ``command.cancel`` frame.

    Parameters
    ----------
    request : int
        The number of the request to cancel.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame(
        {"type": "command.cancel", "request": request},
        version=protocol.CANCEL_PROTOCOL_VERSION,
    )


def commands_list_frame(request: int) -> bytes:
    """Return a ``commands.list`` frame.

    Parameters
    ----------
    request : int
        The request number.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame(
        {"type": "commands.list", "request": request},
        version=protocol.COMMAND_PROTOCOL_VERSION,
    )


def command_frame(request: int, command: str, input: Mapping[str, Any]) -> bytes:
    """Return a ``command.run`` frame.

    The function copies the input mapping. It does not change the mapping of the
    caller. The input must be JSON values. A NaN or an infinite number raises an error.

    Parameters
    ----------
    request : int
        The request number.
    command : str
        The name of the command, for example ``"object.set_visible"``.
    input : Mapping[str, Any]
        The arguments of the command.

    Returns
    -------
    bytes
        The frame.

    Raises
    ------
    ProtocolError
        If the input cannot be encoded as JSON, for example for a NaN value.
    """
    header = {
        "type": "command.run",
        "request": request,
        "command": command,
        "input": dict(input),
    }
    try:
        return protocol.encode_frame(header, version=protocol.COMMAND_PROTOCOL_VERSION)
    except ValueError as error:
        raise protocol.ProtocolError(
            "The input of the command has a value that JSON cannot store, such as NaN "
            "or an infinite number. Use finite numbers."
        ) from error


def _upsert_header(
    id: str,
    object_type: str,
    name: str | None,
    parent: str | None,
    revision: int,
) -> dict[str, Any]:
    """Return the header of a ``scene.upsert`` frame without data."""
    return {
        "type": "scene.upsert",
        "id": id,
        "object_type": object_type,
        "name": id if name is None else name,
        "parent": parent,
        "revision": revision,
        "edit": [],
        "replaces": None,
        "geometry": {},
        "params": {},
        "data": [],
        "property_groups": [],
        "drape": None,
        "style": {},
    }


def group_frame(
    id: str,
    name: str | None,
    parent: str | None,
    revision: int,
) -> bytes:
    """Return a ``scene.upsert`` frame for a group.

    Parameters
    ----------
    id : str
        The identifier of the group.
    name : str or None
        The display name. If it is None, the name is the identifier.
    parent : str or None
        The identifier of the parent group.
    revision : int
        The revision number.

    Returns
    -------
    bytes
        The frame.
    """
    return protocol.encode_frame(_upsert_header(id, "Group", name, parent, revision))


def _utf8_size(value_map: list[dict[str, Any]], name: str) -> int:
    """Return the declared bytes of a value map. Raise if a name is not valid UTF-8."""
    total = 0
    try:
        for item in value_map:
            total += 8 + len(item["name"].encode("utf-8"))
    except UnicodeEncodeError as error:
        raise InvalidInputError(
            f"A category name of the property '{name}' cannot be encoded as UTF-8. "
            "Remove the characters that are not valid."
        ) from error
    return total


def _declared_bytes(name: str, value: object) -> int:
    """Return the bytes that a property needs, from its shape and target type.

    The function does not copy the values. For text it returns an upper bound.
    """
    if isinstance(value, Categorical):
        codes = np.asarray(value.codes)
        total = int(codes.size) * 4
        for label in value.names:
            try:
                total += 8 + len(str(label).encode("utf-8"))
            except UnicodeEncodeError as error:
                raise InvalidInputError(
                    f"A category name of the property '{name}' cannot be encoded as UTF-8. "
                    "Remove the characters that are not valid."
                ) from error
        return total
    arr = np.asarray(value)
    kind = arr.dtype.kind
    if kind == "f":
        return int(arr.size) * (4 if arr.dtype == np.dtype("float32") else 8)
    if kind == "b":
        return int(arr.size)
    if kind in "iu":
        return int(arr.size) * 4
    if kind in "US":
        return int(arr.size) * int(arr.dtype.itemsize) + (int(arr.size) + 1) * 4
    return 0


def _data_array(
    name: str, value: object, association: str, log_set: str | None = None
) -> tuple[dict[str, Any], Buffers]:
    """Return the data description and the buffers of one property.

    Parameters
    ----------
    name : str
        The name of the property.
    value : numpy.ndarray or Categorical
        The values.
    association : str
        ``"vertex"``, ``"cell"`` or ``"interval"``. For a Drillhole, ``"vertex"`` has one
        value for each survey station and ``"interval"`` has one value for each row of
        ``intervals``. An ``"interval"`` array belongs to a log set. ``"cell"`` is for
        the other object types.
    log_set : str or None
        The name of the log set of an ``"interval"`` array. The function writes it to
        the data description as ``log_set``. Use ``None`` for the other associations.

    Returns
    -------
    desc : dict
        The entry for the ``data`` list of the header.
    buffers : list of (str, numpy.ndarray)
        The buffers of the property.

    Raises
    ------
    InvalidInputError
        If the type of the values is not supported, or an integer is outside i32.
    """
    buffer_name = f"data/{name}"
    desc: dict[str, Any] = {
        "id": name,
        "name": name,
        "buffer": buffer_name,
        "association": association,
    }
    if log_set is not None:
        desc["log_set"] = log_set
    if isinstance(value, Categorical):
        codes = np.asarray(value.codes).astype(np.uint32, copy=False)
        colours = value.colours
        value_map = []
        for index, label in enumerate(value.names):
            if colours is None:
                colour = list(DEFAULT_CATEGORY_COLOURS[index % len(DEFAULT_CATEGORY_COLOURS)])
            else:
                colour = [int(c) for c in colours[index]]
            value_map.append({"id": index, "name": str(label), "colour": colour})
        _utf8_size(value_map, name)
        desc["data_type"] = "referenced"
        desc["value_map"] = value_map
        return desc, [(buffer_name, codes)]
    arr = np.asarray(value)
    kind = arr.dtype.kind
    if kind == "f":
        if arr.dtype != np.dtype("float32"):
            arr = arr.astype(np.float64, copy=False)
        desc["data_type"] = "float"
        return desc, [(buffer_name, arr)]
    if kind == "b":
        desc["data_type"] = "boolean"
        return desc, [(buffer_name, arr.astype(np.uint8, copy=False))]
    if kind in "iu":
        if arr.size and (arr.min() < _I32_MIN or arr.max() > _I32_MAX):
            raise InvalidInputError(
                f"The property '{name}' has integers outside the i32 range. "
                "Use values from -2147483648 to 2147483647."
            )
        desc["data_type"] = "integer"
        return desc, [(buffer_name, arr.astype(np.int32, copy=False))]
    if kind in "US":
        try:
            if kind == "U":
                parts = [str(item).encode("utf-8") for item in arr.tolist()]
            else:
                parts = [bytes(item) for item in arr.tolist()]
                for part in parts:
                    part.decode("utf-8")
        except (UnicodeEncodeError, UnicodeDecodeError) as error:
            raise InvalidInputError(
                f"The property '{name}' has text that is not valid UTF-8. "
                "Remove the characters that are not valid."
            ) from error
        offsets = np.zeros(len(parts) + 1, dtype=np.uint32)
        if parts:
            offsets[1:] = np.cumsum([len(p) for p in parts], dtype=np.uint64)
        data = np.frombuffer(b"".join(parts), dtype=np.uint8)
        offsets_name = f"{buffer_name}/offsets"
        desc["data_type"] = "text"
        desc["offsets"] = offsets_name
        return desc, [(buffer_name, data), (offsets_name, offsets)]
    raise InvalidInputError(
        f"The property '{name}' has the type {arr.dtype}, which is not supported. "
        "Use numbers, booleans, text or Categorical."
    )


class _Builder:
    """Collect the header parts and the buffers of one object."""

    def __init__(self, id: str, origin: np.ndarray) -> None:
        self.id = id
        self.origin = origin
        self.buffers: Buffers = []
        self.names: set[str] = set()
        self.geometry: dict[str, str] = {}
        self.data: list[dict[str, Any]] = []
        self.groups: list[dict[str, Any]] = []
        self.declared = 0

    def reserve(self, size: int) -> None:
        """Add a declared size. Raise if a limit is passed. Call before any copy."""
        self.declared += size
        limit = protocol.MAX_OBJECT_BYTES
        if self.declared > limit:
            raise ObjectTooLargeError(self.id, self.declared, limit)
        frame_limit = protocol.MAX_FRAME_BYTES
        if self.declared > frame_limit - protocol.MAX_HEADER_BYTES:
            raise ObjectTooLargeError(self.id, self.declared, frame_limit)

    def add_buffer(self, name: str, array: np.ndarray) -> None:
        """Add a buffer. Raise if the name is already used."""
        if name in self.names:
            raise InvalidInputError(
                f"The buffer name '{name}' is used two times. Use different property names."
            )
        self.names.add(name)
        self.buffers.append((name, array))

    def add_geometry(self, role: str, array: np.ndarray) -> None:
        """Add a geometry buffer for a role."""
        name = f"geometry/{role}"
        self.add_buffer(name, array)
        self.geometry[role] = name

    def add_vertices(self, verts: np.ndarray) -> None:
        """Shift the vertices in f64 and add them as f32 or f64."""
        # Use the declared size first. Do not copy the array before this check.
        self.reserve(int(verts.size) * 8)
        # The shift makes one new f64 array. The array of the caller does not change.
        shifted = np.subtract(verts, self.origin, dtype=np.float64)
        low = float(shifted.min()) if shifted.size else 0.0
        high = float(shifted.max()) if shifted.size else 0.0
        if not (np.isfinite(low) and np.isfinite(high)):
            raise InvalidInputError("The vertices must be finite numbers. Check the input.")
        if max(abs(low), abs(high)) > protocol.F32_LIMIT:
            self.add_geometry("vertices", shifted)
        else:
            self.add_geometry("vertices", shifted.astype(np.float32))

    def add_index(self, role: str, array: object) -> None:
        """Add an index buffer as u32."""
        view = np.asarray(array)
        self.reserve(int(view.size) * 4)
        self.add_geometry(role, view.astype(np.uint32, copy=False))

    def add_property(
        self, name: str, value: object, association: str, log_set: str | None = None
    ) -> None:
        """Add a property with its buffers. A log set name is optional."""
        self.reserve(_declared_bytes(name, value))
        desc, buffers = _data_array(name, value, association, log_set)
        for buffer_name, array in buffers:
            self.add_buffer(buffer_name, array)
        self.data.append(desc)

    def add_properties(
        self, props: Mapping[str, object], association: str, log_set: str | None = None
    ) -> None:
        """Add a mapping of properties. A log set name is optional."""
        for name, value in props.items():
            self.add_property(name, value, association, log_set)

    def add_pair_group(self, role: str, names: tuple[str, str], pair: Sequence[object]) -> None:
        """Add a property group of two float arrays, such as orientation."""
        members = []
        for member_role, values in zip(names, pair):
            data_name = f"{role}/{member_role}"
            self.add_property(data_name, np.asarray(values, dtype=np.float64), "vertex")
            members.append({"role": member_role, "data": data_name})
        self.groups.append({"id": role, "name": role, "role": role, "members": members})


def _finite_f64(value: object, label: str) -> np.ndarray:
    """Return a finite float64 array. Raise if a value is NaN or infinite."""
    arr = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise InvalidInputError(f"The field '{label}' must be finite numbers. Check the input.")
    return arr


def encode_object(
    obj: SceneObject,
    *,
    id: str,
    name: str | None,
    parent: str | None,
    revision: int,
    origin: Sequence[float],
) -> bytes:
    """Return the ``scene.upsert`` frame for an object.

    The function subtracts the origin from the coordinates in float64. Then it
    writes the vertices as f32 if all values are near the origin, or as f64 if not.

    Parameters
    ----------
    obj : SceneObject
        The object to encode.
    id : str
        The identifier of the object.
    name : str or None
        The display name. If it is None, the name is the identifier.
    parent : str or None
        The identifier of the parent group.
    revision : int
        The revision number.
    origin : sequence of float
        The scene origin, three values.

    Returns
    -------
    bytes
        The frame.

    Raises
    ------
    InvalidInputError
        If a coordinate is NaN or infinite, or a property is not supported.
    ObjectTooLargeError
        If the object or the frame is over the size limit.
    """
    try:
        origin_arr = np.asarray(origin, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise InvalidInputError("The origin must be three numbers. Check the input.") from error
    if origin_arr.shape != (3,) or not np.all(np.isfinite(origin_arr)):
        raise InvalidInputError("The origin must be three finite numbers. Check the input.")
    builder = _Builder(id, origin_arr)
    params: dict[str, Any] = {}
    log_sets: list[dict[str, str]] = []
    if isinstance(obj, Drillhole):
        object_type = "Drillhole"
        builder.reserve(
            (3 + int(np.asarray(obj.survey).size) + int(np.asarray(obj.intervals).size)) * 8
        )
        collar = _finite_f64(obj.collar, "collar") - origin_arr
        builder.add_geometry("collar", collar)
        builder.add_geometry("survey", _finite_f64(obj.survey, "survey"))
        builder.add_geometry("intervals", _finite_f64(obj.intervals, "intervals"))
        if obj.end_of_hole is not None:
            params["end_of_hole"] = [float(obj.end_of_hole)]
        builder.add_properties(obj.properties, "vertex")
        if obj.interval_properties:
            log_sets.append({"name": _INTERVAL_LOG_SET, "kind": "interval"})
            builder.add_properties(obj.interval_properties, "interval", _INTERVAL_LOG_SET)
    elif isinstance(obj, BlockModel):
        object_type = "BlockModel"
        shifted = _finite_f64(obj.origin, "origin") - origin_arr
        params["origin"] = [float(v) for v in shifted]
        params["cell_sizes"] = [float(v) for v in np.asarray(obj.cell_sizes)]
        params["cells"] = [int(v) for v in np.asarray(obj.cells)]
        builder.add_properties(obj.properties, "vertex")
        builder.add_properties(obj.cell_properties, "cell")
    else:
        if isinstance(obj, Lines):
            object_type, index_role, index_attr = "Curve", "segments", "segments"
        elif isinstance(obj, Surface):
            object_type, index_role, index_attr = "Surface", "triangles", "triangles"
        elif isinstance(obj, TetraMesh):
            object_type, index_role, index_attr = "TetraMesh", "tetrahedra", "tetrahedra"
        elif isinstance(obj, (Points, Vectors)):
            object_type, index_role, index_attr = "Points", "", ""
        else:
            raise InvalidInputError(
                f"The object type {type(obj).__name__} is not supported. "
                "Use Points, Vectors, Lines, Surface, BlockModel, TetraMesh or Drillhole."
            )
        builder.add_vertices(obj.vertices)
        if index_role:
            builder.add_index(index_role, getattr(obj, index_attr))
        if isinstance(obj, Vectors):
            if obj.orientation is not None:
                builder.add_pair_group("orientation", ("dip_direction", "dip"), obj.orientation)
            if obj.lineation is not None:
                builder.add_pair_group("lineation", ("trend", "plunge"), obj.lineation)
        builder.add_properties(obj.properties, "vertex")
        cell_properties = getattr(obj, "cell_properties", None)
        if cell_properties:
            builder.add_properties(cell_properties, "cell")

    header = _upsert_header(id, object_type, name, parent, revision)
    header["geometry"] = builder.geometry
    header["params"] = params
    header["data"] = builder.data
    header["property_groups"] = builder.groups
    if log_sets:
        header["log_sets"] = log_sets
    try:
        return protocol.encode_frame(header, builder.buffers)
    except protocol.ProtocolError as error:
        raise ObjectTooLargeError(id, builder.declared, protocol.MAX_FRAME_BYTES) from error
