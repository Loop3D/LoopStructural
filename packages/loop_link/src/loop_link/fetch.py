"""Fetch viewer objects and decode them into scene objects.

The viewer sends an object as a JSON header and a set of named arrays. The
function :func:`decode_object` makes a scene object from them. The coordinates
of the result are world coordinates in float64. The function
:func:`fetch_object` asks the viewer for an object and collects the reply.
"""

from __future__ import annotations

import itertools
import math
import threading
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from loop_link import encode, protocol
from loop_link.errors import (
    InvalidInputError,
    ObjectTooLargeError,
    UnsupportedObjectError,
    ViewerClosedError,
)
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
from loop_link.protocol import Frame, ProtocolError
from loop_link.tool import ToolConnection

__all__ = ["FetchResult", "decode_object", "fetch_object"]

# The numpy type of each array type name in the protocol. All are little-endian.
_ARRAY_DTYPES = {
    "f32": np.dtype("<f4"),
    "f64": np.dtype("<f8"),
    "u32": np.dtype("<u4"),
    "i32": np.dtype("<i4"),
    "u8": np.dtype("u1"),
}

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


def _is_int(value: object) -> bool:
    """Return whether the value is an integer and not a boolean."""
    return isinstance(value, int) and not isinstance(value, bool)


def _as_dict(value: object, what: str) -> dict[str, Any]:
    """Return the value if it is a dict, or raise a ProtocolError."""
    if not isinstance(value, dict):
        raise ProtocolError(f"The {what} is not an object.")
    return value


def _as_origin(origin: Sequence[float]) -> tuple[float, float, float]:
    """Return the origin as a tuple of 3 floats."""
    x, y, z = (float(value) for value in origin)
    return (x, y, z)


class _Assembler:
    """Collect the frames of one ``object.fetch`` reply and make the result.

    Parameters
    ----------
    object_id : str
        The identifier of the object that was requested.
    """

    def __init__(self, object_id: str) -> None:
        self._object_id = object_id
        self._header: dict[str, Any] | None = None
        self._arrays: dict[str, np.ndarray] = {}
        self._shapes: dict[str, tuple[int, ...]] = {}
        self._dtypes: dict[str, np.dtype] = {}
        self._filled: dict[str, int] = {}
        self._parts = 0
        self._ended = False

    def feed(self, frame: Frame) -> None:
        """Take the next frame of the reply.

        Parameters
        ----------
        frame : loop_link.protocol.Frame
            An ``object.header``, ``object.part`` or ``object.end`` frame.

        Raises
        ------
        loop_link.protocol.ProtocolError
            If the frame breaks a rule of the stream.
        loop_link.errors.ObjectTooLargeError
            If the arrays of the header are larger than the limit.
        """
        kind = frame.header.get("type")
        if self._ended:
            raise ProtocolError(f"The viewer sent a '{kind}' frame after 'object.end'.")
        if self._header is None:
            if kind != "object.header":
                raise ProtocolError(
                    f"The first frame of an object must be 'object.header', not '{kind}'."
                )
            self._take_header(frame.header)
        elif kind == "object.part":
            self._take_part(frame)
        elif kind == "object.end":
            self._take_end(frame.header)
        else:
            raise ProtocolError(f"The viewer sent an unknown frame '{kind}' for an object.")

    def finish(self) -> FetchResult:
        """Return the object that the frames described.

        Returns
        -------
        FetchResult
            The decoded object and its data from the header.

        Raises
        ------
        loop_link.protocol.ProtocolError
            If the stream did not end, or the header is not valid.
        """
        header = self._header
        if header is None or not self._ended:
            raise ProtocolError("The viewer did not finish the object.")
        arrays = {
            name: data.view(self._dtypes[name]).reshape(self._shapes[name])
            for name, data in self._arrays.items()
        }
        try:
            obj = header["object"]
            origin = header["origin"]
            scene_object = decode_object(obj, arrays, origin)
            return FetchResult(
                id=self._object_id,
                object=scene_object,
                object_type=str(obj["object_type"]),
                revision=int(obj["revision"]),
                name=obj.get("name"),
                parent=obj.get("parent"),
                origin=None if origin is None else _as_origin(origin),
                crs=header.get("crs"),
            )
        except (KeyError, TypeError, ValueError, IndexError) as error:
            raise ProtocolError(f"The object header is not valid: {error!r}") from error

    def _take_header(self, header: dict[str, Any]) -> None:
        """Check the header, then make the empty arrays."""
        entries = header.get("arrays")
        if not isinstance(entries, list):
            raise ProtocolError("The 'arrays' field of the object header is not a list.")
        sizes: dict[str, int] = {}
        shapes: dict[str, tuple[int, ...]] = {}
        dtypes: dict[str, np.dtype] = {}
        for entry in entries:
            entry = _as_dict(entry, "array entry")
            name, dtype_name, shape = entry.get("name"), entry.get("dtype"), entry.get("shape")
            if not isinstance(name, str) or name in sizes:
                raise ProtocolError(f"The array name {name!r} is not valid or not unique.")
            if not isinstance(dtype_name, str) or dtype_name not in _ARRAY_DTYPES:
                raise ProtocolError(f"The array '{name}' has the type {dtype_name!r}.")
            if not isinstance(shape, list) or not all(_is_int(n) and n >= 0 for n in shape):
                raise ProtocolError(f"The array '{name}' has a shape that is not valid.")
            dtypes[name] = _ARRAY_DTYPES[dtype_name]
            shapes[name] = tuple(shape)
            sizes[name] = math.prod(shape) * dtypes[name].itemsize
        total = sum(sizes.values())
        if total > protocol.MAX_OBJECT_BYTES:
            raise ObjectTooLargeError(self._object_id, total, protocol.MAX_OBJECT_BYTES)
        obj = _as_dict(header.get("object"), "object header")
        self._check_buffer_names(obj, sizes)
        self._header = header
        self._shapes = shapes
        self._dtypes = dtypes
        self._arrays = {name: np.empty(size, dtype=np.uint8) for name, size in sizes.items()}
        self._filled = dict.fromkeys(sizes, 0)

    @staticmethod
    def _check_buffer_names(obj: dict[str, Any], sizes: Mapping[str, int]) -> None:
        """Check that each buffer name of the object is in the arrays."""
        names = list(_as_dict(obj.get("geometry", {}), "geometry").values())
        data = obj.get("data", [])
        if not isinstance(data, list):
            raise ProtocolError("The 'data' field of the object is not a list.")
        for desc in data:
            desc = _as_dict(desc, "data entry")
            names.append(desc.get("buffer"))
            if "offsets" in desc:
                names.append(desc["offsets"])
        for name in names:
            if not isinstance(name, str) or name not in sizes:
                raise ProtocolError(f"The object uses the buffer {name!r}, which is not sent.")

    def _take_part(self, frame: Frame) -> None:
        """Copy the bytes of one ``object.part`` frame into its array."""
        header = frame.header
        name, offset, size = header.get("array"), header.get("offset"), header.get("bytes")
        if not isinstance(name, str) or name not in self._arrays:
            raise ProtocolError(f"The viewer sent a part for the unknown array {name!r}.")
        if not _is_int(offset) or not _is_int(size) or offset < 0 or size < 0:
            raise ProtocolError(f"The part of the array '{name}' has a bad offset or size.")
        if offset % 8 != 0:
            raise ProtocolError(
                f"The part of the array '{name}' starts at {offset}, not a multiple of 8."
            )
        part = frame.buffers.get("part")
        if not isinstance(part, np.ndarray) or part.dtype != np.uint8 or part.size != size:
            raise ProtocolError(f"The part of the array '{name}' is not {size} bytes of type u8.")
        if offset != self._filled[name]:
            raise ProtocolError(
                f"The part of the array '{name}' starts at {offset}, "
                f"but {self._filled[name]} bytes are filled."
            )
        target = self._arrays[name]
        if offset + size > target.size:
            raise ProtocolError(f"The part of the array '{name}' goes past the end of the array.")
        target[offset : offset + size] = part.reshape(-1)
        self._filled[name] = offset + size
        self._parts += 1

    def _take_end(self, header: dict[str, Any]) -> None:
        """Check the ``object.end`` frame."""
        parts = header.get("parts")
        if not _is_int(parts) or parts != self._parts:
            raise ProtocolError(f"The viewer said {parts!r} parts, but sent {self._parts}.")
        for name, array in self._arrays.items():
            if self._filled[name] != array.size:
                raise ProtocolError(f"The array '{name}' is not full when the object ends.")
        self._ended = True


def fetch_object(
    tool: ToolConnection,
    id: str,
    *,
    progress: Callable[[float], None] | None = None,
    cancel: threading.Event | None = None,
    timeout: float = 30.0,
) -> FetchResult:
    """Ask the viewer for one object and decode it.

    Parameters
    ----------
    tool : loop_link.tool.ToolConnection
        The connection to the viewer.
    id : str
        The identifier of the object.
    progress : callable, optional
        Called with the fraction of each progress report.
    cancel : threading.Event, optional
        When it is set, the request is cancelled.
    timeout : float, optional
        The time in seconds to wait for the next frame.

    Returns
    -------
    FetchResult
        The object, in world coordinates.

    Raises
    ------
    loop_link.errors.ViewerClosedError
        If the connection closed, or the viewer sent an object that is not valid.
        In the second case the function closes the connection.
    loop_link.errors.ObjectTooLargeError
        If the object is larger than the limit.
    """
    assembler = _Assembler(id)
    try:
        tool.request(
            lambda n: encode.fetch_frame(n, id),
            until=("object.end",),
            feature="selection and fetch (geoviewer protocol 1.18)",
            object_id=id,
            on_frame=assembler.feed,
            progress=progress,
            cancel=cancel,
            timeout=timeout,
        )
        return assembler.finish()
    except ProtocolError as error:
        tool.close()
        raise ViewerClosedError(f"The viewer sent an object that is not valid: {error}") from error
