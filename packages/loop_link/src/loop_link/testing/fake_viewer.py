"""A fake viewer for the tests: a WebSocket server on 127.0.0.1.

The fake viewer does what the link of the real viewer does. It checks the bearer
token, answers ``hello`` with ``hello`` and ``scene.inventory``, and keeps the
objects that it receives. It has its own small server-side WebSocket code.
"""

from __future__ import annotations

import base64
import copy
import hashlib
import json
import math
import os
import socket
import struct
import threading
import time
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .. import protocol
from ..protocol import (
    MAX_FRAME_BYTES,
    MAX_HEADER_BYTES,
    MAX_OBJECT_BYTES,
    Frame,
    ProtocolError,
)

_GUID = b"258EAFA5-E914-47DA-95CA-C5AB0DC85B11"

# The geometry roles that each object type needs.
_REQUIRED_ROLES: dict[str, tuple[str, ...]] = {
    "Points": ("vertices",),
    "Curve": ("vertices", "segments"),
    "Surface": ("vertices", "triangles"),
    "TetraMesh": ("vertices", "tetrahedra"),
    "BlockModel": (),
    "Drillhole": ("collar", "survey", "intervals"),
    "Group": (),
}
# The roles whose rows hold vertex indices, and the roles that hold coordinates.
_INDEX_ROLES = ("segments", "triangles", "tetrahedra")
_COORDINATE_ROLES = ("vertices", "collar", "survey")
# The roles that give the cell count of an object.
_CELL_ROLES = ("segments", "triangles", "tetrahedra")
# The row width that each geometry role needs.
_ROLE_WIDTH = {
    "vertices": 3,
    "survey": 3,
    "intervals": 2,
    "segments": 2,
    "triangles": 3,
    "tetrahedra": 4,
}


# The protocol names of the array types, by NumPy dtype name.
_PROTOCOL_DTYPES = {
    "float32": "f32",
    "float64": "f64",
    "uint32": "u32",
    "int32": "i32",
    "uint8": "u8",
}


class _Closed(Exception):
    """The peer closed the connection."""


@dataclass
class StoredObject:
    """An object that the fake viewer keeps.

    Parameters
    ----------
    header : dict
        The header of the ``scene.upsert`` frame.
    buffers : dict of str to numpy.ndarray
        The arrays of the frame, by name.
    """

    header: dict[str, Any]
    buffers: dict[str, np.ndarray] = field(default_factory=dict)


def _error(code: str, message: str, object_id: object = None) -> Frame:
    """Return an ``error`` frame."""
    header: dict[str, Any] = {"type": "error", "code": code, "message": message}
    if object_id is not None:
        header["object"] = object_id
    return Frame(header=header, buffers={})


def _encode(frame: Frame, version: dict[str, int] | None = None) -> bytes:
    """Return the bytes of a frame, with ``version`` in the header if it is given."""
    return protocol.encode_frame(frame.header, list(frame.buffers.items()), version=version)


def _request_error(
    request: object, code: str, message: str, field_name: str | None = None
) -> Frame:
    """Return a ``request.error`` frame for a tool request."""
    error: dict[str, Any] = {"code": code, "message": message}
    if field_name is not None:
        error["field"] = field_name
    return Frame(header={"type": "request.error", "request": request, "error": error}, buffers={})


# The kinds of selection that each object type has. Drillhole and Group have none.
_SELECTION_KINDS: dict[str, tuple[str, ...]] = {
    "Points": ("nodes",),
    "Curve": ("nodes", "cells"),
    "Surface": ("nodes", "cells"),
    "BlockModel": ("nodes", "cells"),
    "TetraMesh": ("nodes", "cells"),
    "Drillhole": (),
    "Group": (),
}


# The value types of the command inputs, with the JSON schema of each type.
_INPUT_SCHEMAS: dict[str, dict[str, Any]] = {
    "str": {"type": "string"},
    "str_or_null": {"type": ["string", "null"]},
    "bool": {"type": "boolean"},
    "u32_list": {"type": "array", "items": {"type": "integer", "minimum": 0}},
    "colour": {
        "type": ["array", "null"],
        "items": {"type": "integer", "minimum": 0, "maximum": 255},
        "minItems": 3,
        "maxItems": 3,
    },
    "scale": {"type": "string", "enum": ["linear", "log"]},
    "edit": {"type": "object"},
}
# The commands: (description, needs_window, {field: (value type, required)}).
_COMMANDS: dict[str, tuple[str, bool, dict[str, tuple[str, bool]]]] = {
    "filter.edit": (
        "Change the filters of an object.",
        False,
        {"object": ("str", True), "edit": ("edit", True)},
    ),
    "colour.set_active": (
        "Choose the array that colours an object.",
        False,
        {"id": ("str", True), "array": ("str_or_null", False)},
    ),
    "colour.set_map": (
        "Choose the colour map of an object.",
        False,
        {"id": ("str", True), "name": ("str", True)},
    ),
    "colour.set_range": (
        "Set the colour range of an object.",
        False,
        {"id": ("str", True), "min": ("str", True), "max": ("str", True)},
    ),
    "colour.clear_range": (
        "Clear the colour range of an object.",
        False,
        {"id": ("str", True)},
    ),
    "colour.set_scale": (
        "Set the colour scale of an object.",
        True,
        {"object": ("str", True), "scale": ("scale", True)},
    ),
    "object.set_colour": (
        "Set the plain colour of an object, or clear it.",
        False,
        {"id": ("str", True), "colour": ("colour", True)},
    ),
}
# The fields of each kind of filter edit.
_EDIT_FIELDS: dict[str, dict[str, str]] = {
    "add": {"array": "str"},
    "set_range": {"array": "str", "min": "str", "max": "str", "keep_nan": "bool"},
    "set_values": {"array": "str", "values": "u32_list"},
    "set_on": {"array": "str", "on": "bool"},
    "set_exclude": {"array": "str", "exclude": "bool"},
    "remove": {"array": "str"},
    "cancel": {},
}
# The colour maps that the viewer lists from minor version 18.
_COLOUR_MAPS = ["viridis", "grey", "coolwarm"]
# The data types that have no colour map, range or scale.
_FIXED_COLOUR_TYPES = ("referenced", "boolean")


def _value_ok(kind: str, value: object) -> bool:
    """Return whether ``value`` fits the input type ``kind``."""
    if kind == "str":
        return isinstance(value, str)
    if kind == "str_or_null":
        return value is None or isinstance(value, str)
    if kind == "bool":
        return isinstance(value, bool)
    if kind == "u32_list":
        return isinstance(value, list) and all(type(v) is int and 0 <= v < 1 << 32 for v in value)
    if kind == "colour":
        return value is None or (
            isinstance(value, list)
            and len(value) == 3
            and all(type(v) is int and 0 <= v <= 255 for v in value)
        )
    if kind == "scale":
        return value in ("linear", "log")
    return isinstance(value, dict)


def _check_input(fields: dict[str, tuple[str, bool]], given: object) -> str | None:
    """Return the name of the first bad field of a command input, or ``None``.

    A field is bad if it is unknown, missing, or of the wrong type.
    """
    if not isinstance(given, dict):
        return "input"
    for key in given:
        if key not in fields:
            return str(key)
    for key, (kind, required) in fields.items():
        if key not in given:
            if required:
                return key
        elif not _value_ok(kind, given[key]):
            return key
    return None


def _check_edit(edit: dict[str, Any]) -> str | None:
    """Return the name of the first bad field of a filter edit, or ``None``."""
    kind = edit.get("kind")
    if not isinstance(kind, str) or kind not in _EDIT_FIELDS:
        return "edit"
    fields = {key: (value_kind, True) for key, value_kind in _EDIT_FIELDS[kind].items()}
    fields["kind"] = ("str", True)
    return _check_input(fields, edit)


def _send(conn: socket.socket, payload: bytes) -> None:
    """Write one unmasked binary WebSocket frame."""
    length = len(payload)
    if length < 126:
        head = bytes([0x82, length])
    elif length < 1 << 16:
        head = bytes([0x82, 126]) + struct.pack(">H", length)
    else:
        head = bytes([0x82, 127]) + struct.pack(">Q", length)
    conn.sendall(head + payload)


class _Reader:
    """Read masked client frames from a socket.

    Attributes
    ----------
    deferred : list of bytes
        Messages that a read-ahead took from the socket but that the viewer must
        still handle, in arrival order.
    """

    def __init__(self, conn: socket.socket, buffer: bytes = b"") -> None:
        self._conn = conn
        self._buffer = buffer
        self.deferred: list[bytes] = []

    def _take(self, count: int) -> bytes:
        while len(self._buffer) < count:
            chunk = self._conn.recv(65536)
            if not chunk:
                raise _Closed
            self._buffer += chunk
        data, self._buffer = self._buffer[:count], self._buffer[count:]
        return data

    def _frame(self) -> tuple[bool, int, bytes]:
        """Return (fin, opcode, payload) of one WebSocket frame."""
        head = self._take(2)
        fin = bool(head[0] & 0x80)
        opcode = head[0] & 0x0F
        length = head[1] & 0x7F
        if not head[1] & 0x80:
            raise _Closed
        if length == 126:
            (length,) = struct.unpack(">H", self._take(2))
        elif length == 127:
            (length,) = struct.unpack(">Q", self._take(8))
        if length > MAX_FRAME_BYTES + MAX_HEADER_BYTES:
            raise _Closed
        key = self._take(4)
        payload = self._take(length)
        if length:
            mask = np.frombuffer((key * (length // 4 + 1))[:length], dtype=np.uint8)
            payload = (np.frombuffer(payload, dtype=np.uint8) ^ mask).tobytes()
        return fin, opcode, payload

    def message(self) -> bytes | None:
        """Return the next binary message, or ``None`` after a close frame.

        The method joins fragments and answers ping frames. It ignores text frames.
        """
        parts: list[bytes] = []
        binary = False
        while True:
            fin, opcode, payload = self._frame()
            if opcode == 0x8:
                self._conn.sendall(bytes([0x88, 0]))
                return None
            if opcode == 0x9:
                self._conn.sendall(bytes([0x8A, len(payload[:125])]) + payload[:125])
                continue
            if opcode == 0xA:
                continue
            if opcode in (0x1, 0x2):
                parts = [payload]
                binary = opcode == 0x2
            elif opcode == 0x0:
                parts.append(payload)
            if fin:
                if binary:
                    return b"".join(parts)
                parts = []


def _cell_count(arrays: dict[str, np.ndarray]) -> int | None:
    """Return the number of segments, triangles or tetrahedra, or ``None``."""
    return next((arrays[r].shape[0] for r in _CELL_ROLES if r in arrays), None)


def _expected_length(
    object_type: str,
    association: object,
    arrays: dict[str, np.ndarray],
    params: object,
    log_sets: object,
) -> int | None:
    """Return the length that a data array needs, or ``None`` if it cannot fit."""
    if object_type == "BlockModel":
        cells = params.get("cells") if isinstance(params, dict) else None
        if not isinstance(cells, list) or not all(type(n) is int and n >= 0 for n in cells):
            return None
        if association == "cell":
            return math.prod(cells)
        if association == "vertex":
            return math.prod(n + 1 for n in cells)
        return None
    if object_type == "Drillhole":
        if association == "vertex" and "survey" in arrays:
            return arrays["survey"].shape[0]
        has_set = isinstance(log_sets, list) and any(
            isinstance(item, dict) and item.get("kind") == "interval" for item in log_sets
        )
        if association == "interval" and has_set and "intervals" in arrays:
            return arrays["intervals"].shape[0]
        return None
    vertices = arrays.get("vertices")
    if association == "vertex" and vertices is not None:
        return vertices.shape[0]
    if association == "cell" and object_type in ("Curve", "Surface", "TetraMesh"):
        return _cell_count(arrays)
    return None


def _check_data(
    object_type: str,
    entry: object,
    arrays: dict[str, np.ndarray],
    buffers: dict[str, np.ndarray],
    params: object,
    log_sets: object = None,
) -> str | None:
    """Return the text of a problem in one data array entry, or ``None``.

    An ``interval`` array must name a log set of kind ``interval`` in ``log_set``.
    """
    if not isinstance(entry, dict):
        return "A data entry is not an object."
    name = entry.get("buffer")
    if name not in buffers:
        return f"The data entry {entry.get('name')!r} names the buffer {name!r}, which is missing."
    association = entry.get("association")
    if association not in ("vertex", "cell", "interval"):
        return f"The association {association!r} is not known."
    if association == "interval":
        log_set = entry.get("log_set")
        kinds = {
            item.get("name"): item.get("kind")
            for item in (log_sets if isinstance(log_sets, list) else [])
            if isinstance(item, dict)
        }
        if kinds.get(log_set) != "interval":
            return (
                f"The interval array {name!r} names the log set {log_set!r}. "
                "It needs a log set of kind 'interval' in 'log_sets'."
            )
    length = buffers[name].shape[0] if buffers[name].ndim else 0
    if entry.get("data_type") == "text":
        offsets_name = entry.get("offsets")
        if offsets_name not in buffers:
            return f"The text array {name!r} has no offsets buffer."
        offsets = buffers[offsets_name].astype(np.int64)
        if offsets.ndim != 1 or offsets.size < 1:
            return f"The offsets of {name!r} are not valid."
        if offsets[0] != 0 or offsets[-1] != length or np.any(np.diff(offsets) < 0):
            return f"The offsets of {name!r} do not fit its {length} bytes."
        length = offsets.size - 1
    expected = _expected_length(object_type, association, arrays, params, log_sets)
    if expected is None:
        return f"The association {association!r} does not fit the object type {object_type}."
    if length != expected:
        return f"The data array {name!r} has {length} values. The association needs {expected}."
    return None


_GROUP_MEMBERS: dict[str, tuple[set[str], set[str]]] = {
    "orientation": ({"dip_direction", "dip", "polarity"}, {"dip_direction", "dip"}),
    "lineation": ({"trend", "plunge"}, {"trend", "plunge"}),
}


def _check_group(group: object, data_ids: list[Any]) -> str | None:
    """Return the text of a problem in one property group, or ``None``."""
    if not isinstance(group, dict):
        return "A property group is not an object."
    members = group.get("members") or []
    for member in members:
        if not isinstance(member, dict) or member.get("data") not in data_ids:
            return f"The group {group.get('id')!r} names a data array that the object lacks."
    role = group.get("role")
    if role is None:
        return None
    if role not in _GROUP_MEMBERS:
        return f"The group role {role!r} is not known."
    allowed, required = _GROUP_MEMBERS[role]
    roles = [member.get("role") for member in members]
    if any(r not in allowed for r in roles):
        return (
            f"The {role} group {group.get('id')!r} has a member role that the role does not list."
        )
    if len(set(roles)) != len(roles) or not required <= set(roles):
        return f"The {role} group {group.get('id')!r} needs each of {sorted(required)} once."
    return None


class FakeViewer:
    """A fake viewer that listens on ``127.0.0.1``.

    Parameters
    ----------
    token : str
        The bearer token that a client must send.
    crs : str
        The CRS that the viewer reports in ``hello``.
    origin : tuple of float
        The origin that the viewer reports in ``hello``.
    kind : str
        The kind of the viewer for the instance file, ``"window"`` or other.
    heartbeat_timeout : float
        The seconds without a frame after which a data connection is closed.
    reply_errors : dict of str to str, optional
        Object id to error text. The viewer refuses these objects.
    read_delay : float
        The seconds to wait after each frame.
    protocol_minor : int
        The minor protocol version of the tool endpoint. Below 18, the tool endpoint
        answers each message with an ``unknown_message`` error.
    part_bytes : int
        The size in bytes of one part of an ``object.fetch`` reply.
    fetch_delay : float
        The seconds to wait between the parts of an ``object.fetch`` reply.
    filter_delay : float
        The seconds that the scene stays busy after a ``filter.edit`` that changes a filter.
        For ``protocol_minor`` 20 or more, the seconds that a filter build takes.
    busy_wait : float
        The seconds that a waiting command waits for the busy time to end. After that,
        the viewer answers ``busy``.
    filter_stop_after : float, optional
        For ``protocol_minor`` 20 or more, the seconds after which a filter build gives
        the result ``{"status": "done", "state": "building"}``. The edit stays applied.
        It has no effect if it is not shorter than the build.

    Attributes
    ----------
    port : int
        The port of the server.
    objects : dict of str to StoredObject
        The objects in the scene.
    received : list of Frame
        All frames that the viewer received on all connections.
    tool_requests : list of Frame
        All frames that the viewer received on the tool endpoint.
    """

    def __init__(
        self,
        token: str = "secret",
        crs: str = "EPSG:28350",
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
        kind: str = "window",
        heartbeat_timeout: float = 15.0,
        reply_errors: dict[str, str] | None = None,
        read_delay: float = 0.0,
        protocol_minor: int = 18,
        part_bytes: int = 64 * 1024 * 1024,
        fetch_delay: float = 0.0,
        filter_delay: float = 0.0,
        busy_wait: float = 5.0,
        filter_stop_after: float | None = None,
    ) -> None:
        self.token = token
        self.crs = crs
        self.origin = origin
        self.kind = kind
        self.heartbeat_timeout = heartbeat_timeout
        self.reply_errors = dict(reply_errors or {})
        self.read_delay = read_delay
        self.protocol_minor = protocol_minor
        self.part_bytes = part_bytes
        self.fetch_delay = fetch_delay
        self.filter_delay = filter_delay
        self.busy_wait = busy_wait
        self.filter_stop_after = filter_stop_after
        self._filters: dict[str, dict[str, dict[str, Any]]] = {}
        self._colouring: dict[str, dict[str, Any]] = {}
        self._busy_until = 0.0
        self._selection: dict[str, Any] | None = None
        self.tool_requests: list[Frame] = []
        self.objects: dict[str, StoredObject] = {}
        self.received: list[Frame] = []
        self._lock = threading.Lock()
        self._state_lock = threading.Lock()
        self._data_connections = 0
        self._connections: set[socket.socket] = set()
        self._server = socket.socket()
        self._server.bind(("127.0.0.1", 0))
        self._server.listen(8)
        self.port: int = self._server.getsockname()[1]
        self._thread = threading.Thread(target=self._serve, daemon=True)
        self._thread.start()

    def __enter__(self) -> FakeViewer:
        """Return the viewer."""
        return self

    def __exit__(self, *exc: object) -> None:
        """Close the viewer."""
        self.close()

    def close(self) -> None:
        """Stop the server and close all connections."""
        try:
            self._server.close()
        except OSError:
            pass
        with self._lock:
            connections = list(self._connections)
        for conn in connections:
            try:
                conn.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass
            try:
                conn.close()
            except OSError:
                pass

    # -- state --------------------------------------------------------------

    def inventory(self) -> list[dict[str, Any]]:
        """Return the inventory entries of the objects in the scene."""
        with self._state_lock:
            items = list(self.objects.items())
        return [
            {
                "id": object_id,
                "object_type": stored.header.get("object_type"),
                "revision": stored.header.get("revision", 1),
            }
            for object_id, stored in items
        ]

    @property
    def filters(self) -> dict[str, dict[str, dict[str, Any]]]:
        """A copy of the filters, by object id and then by array name.

        A filter is a dict with the keys ``kind`` (``"range"``, ``"values"`` or ``None``),
        ``on``, ``minimum``, ``maximum``, ``keep_nan``, ``values`` and ``exclude``.
        """
        with self._state_lock:
            return copy.deepcopy(self._filters)

    @property
    def colouring(self) -> dict[str, dict[str, Any]]:
        """A copy of the colouring, by object id.

        Each value is a dict with the keys ``array``, ``colour_map``, ``range``,
        ``scale`` and ``plain``.
        """
        with self._state_lock:
            return copy.deepcopy(self._colouring)

    @property
    def selection(self) -> dict[str, Any] | None:
        """The selection, or ``None``.

        The dict has the keys ``object``, ``revision``, ``cells`` and ``nodes``. The
        value of ``cells`` or ``nodes`` is a ``uint32`` array without repeats, or
        ``None`` if the selection has no index of that kind.
        """
        with self._state_lock:
            return None if self._selection is None else dict(self._selection)

    def set_user_selection(
        self,
        object_id: str,
        *,
        cells: Sequence[int] | None = None,
        nodes: Sequence[int] | None = None,
    ) -> None:
        """Select cells or nodes, as if the user did it in the window.

        Parameters
        ----------
        object_id : str
            The id of the object.
        cells, nodes : sequence of int, optional
            The indices to select.

        Raises
        ------
        ValueError
            If the object is not known, the object lacks a kind, or an index is out of range.
        """
        kinds = {
            kind: np.asarray(values, dtype=np.int64)
            for kind, values in (("cells", cells), ("nodes", nodes))
            if values is not None
        }
        with self._state_lock:
            problem = self._selection_problem(object_id, None, kinds)
            if problem is not None:
                raise ValueError(problem[1])
            self._store_selection(object_id, kinds)

    def _store_selection(self, object_id: str, kinds: dict[str, np.ndarray]) -> None:
        """Replace the selection. The caller holds the state lock."""
        unique = {
            kind: np.unique(values).astype(np.uint32)
            for kind, values in kinds.items()
            if values.size
        }
        if not unique:
            self._selection = None
            return
        self._selection = {
            "object": object_id,
            "revision": self.objects[object_id].header.get("revision", 1),
            "cells": unique.get("cells"),
            "nodes": unique.get("nodes"),
        }

    @staticmethod
    def _element_count(stored: StoredObject, kind: str) -> int:
        """Return the number of nodes or cells of an object."""
        header = stored.header
        if header.get("object_type") == "BlockModel":
            params = header.get("params")
            shape = params.get("cells") if isinstance(params, dict) else None
            shape = [int(n) for n in shape] if isinstance(shape, list) else [0, 0, 0]
            if kind == "cells":
                return math.prod(shape)
            return math.prod(n + 1 for n in shape)
        geometry = header.get("geometry") or {}
        roles = ("vertices",) if kind == "nodes" else _CELL_ROLES
        for role in roles:
            if role in geometry and geometry[role] in stored.buffers:
                return int(stored.buffers[geometry[role]].shape[0])
        return 0

    def _selection_problem(
        self, object_id: object, revision: object, kinds: dict[str, np.ndarray]
    ) -> tuple[str, str] | None:
        """Return (field, message) for the first problem of a selection, or ``None``.

        The caller holds the state lock.
        """
        stored = self.objects.get(object_id) if isinstance(object_id, str) else None
        if stored is None:
            return "object", f"The object {object_id!r} is not in the scene."
        current = stored.header.get("revision", 1)
        if revision is not None and revision != current:
            return "revision", f"The revision {revision!r} is not the revision {current}."
        allowed = _SELECTION_KINDS.get(stored.header.get("object_type"), ())
        for kind, values in kinds.items():
            if kind not in allowed:
                return "kind", f"The object {object_id!r} has no {kind}."
            count = self._element_count(stored, kind)
            if values.size and (int(values.min()) < 0 or int(values.max()) >= count):
                return "indices", f"The {kind} need indices from 0 to {count - 1}."
        return None

    def write_instance_file(self, settings: Path) -> Path:
        """Write ``instances/<pid>.json`` in a settings folder and return its path.

        The method does not write ``link.json``.
        """
        folder = Path(settings) / "instances"
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / f"{os.getpid()}.json"
        data = {"port": self.port, "kind": self.kind, "pid": os.getpid(), "focused": time.time()}
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    # -- frames -------------------------------------------------------------

    def accept_frame(self, data: bytes) -> list[Frame]:
        """Apply one data-client frame and return the reply frames.

        Parameters
        ----------
        data : bytes
            The frame from the client.

        Returns
        -------
        list of Frame
            The frames that the viewer sends back. The list can be empty.
        """
        try:
            frame = protocol.decode_frame(data)
        except ProtocolError as error:
            return [_error("invalid_frame", str(error))]
        with self._state_lock:
            self.received.append(frame)
        header = frame.header
        version = header.get("version")
        major = version.get("major") if isinstance(version, dict) else None
        if major != protocol.PROTOCOL_VERSION["major"]:
            return [_error("unsupported_version", f"The version {version!r} is not supported.")]
        kind = header.get("type")
        if kind == "hello":
            hello = {
                "type": "hello",
                "application": {"name": "fake-viewer"},
                "crs": self.crs,
                "origin": list(self.origin),
            }
            inventory = {"type": "scene.inventory", "entries": self.inventory()}
            return [Frame(header=hello, buffers={}), Frame(header=inventory, buffers={})]
        if kind == "heartbeat":
            return []
        if kind == "scene.upsert":
            return self._upsert(frame)
        if kind == "scene.remove":
            with self._state_lock:
                self.objects.pop(header.get("id"), None)
            return []
        if kind == "scene.clear":
            with self._state_lock:
                self.objects.clear()
            return []
        return [_error("unknown_message", f"The message type {kind!r} is not known.")]

    def _upsert(self, frame: Frame) -> list[Frame]:
        """Check and keep one object."""
        header = frame.header
        object_id = header.get("id")
        if object_id in self.reply_errors:
            return [_error("invalid_object", self.reply_errors[object_id], object_id)]
        problem = self._check_object(frame)
        if problem is not None:
            return [_error("invalid_object", problem, object_id)]
        names = {entry.get("name") for entry in header.get("data") or [] if isinstance(entry, dict)}
        with self._state_lock:
            self.objects[object_id] = StoredObject(header=header, buffers=dict(frame.buffers))
            kept = self._filters.get(object_id)
            if kept is not None:
                self._filters[object_id] = {k: v for k, v in kept.items() if k in names}
        return []

    @staticmethod
    def _check_object(frame: Frame) -> str | None:
        """Return the text of the first problem of an object, or ``None``."""
        header, buffers = frame.header, frame.buffers
        object_id = header.get("id")
        object_type = header.get("object_type")
        if not isinstance(object_id, str):
            return "The object has no id."
        if object_type not in _REQUIRED_ROLES:
            return f"The object type {object_type!r} is not known."
        geometry = header.get("geometry") or {}
        if not isinstance(geometry, dict):
            return "The field 'geometry' is not an object."
        for role in _REQUIRED_ROLES[object_type]:
            if role not in geometry:
                return f"The geometry role {role!r} is missing."
        arrays: dict[str, np.ndarray] = {}
        for role, name in geometry.items():
            if name not in buffers:
                return f"The geometry role {role!r} names the buffer {name!r}, which is missing."
            arrays[role] = buffers[name]
        collar = arrays.get("collar")
        if collar is not None and collar.shape not in ((3,), (1, 3)):
            return f"The buffer for 'collar' has the shape {collar.shape}. It needs (3,) or (1, 3)."
        for role, width in _ROLE_WIDTH.items():
            array = arrays.get(role)
            if array is not None and (array.ndim != 2 or array.shape[1] != width):
                return (
                    f"The buffer for {role!r} has the shape {array.shape}. It needs (n, {width})."
                )
        for role in _COORDINATE_ROLES:
            array = arrays.get(role)
            if array is not None and array.dtype.kind == "f" and not np.all(np.isfinite(array)):
                return f"The buffer for {role!r} has a coordinate that is not finite."
        if "vertices" in arrays:
            vertices = arrays["vertices"]
            count = vertices.shape[0] if vertices.ndim else 0
            for role in _INDEX_ROLES:
                array = arrays.get(role)
                if array is None or not array.size:
                    continue
                if array.dtype.kind not in "iu":
                    return f"The buffer for {role!r} does not hold integers."
                if int(array.min()) < 0 or int(array.max()) >= count:
                    return f"The buffer for {role!r} has an index outside 0 to {count - 1}."
        revision = header.get("revision")
        if type(revision) is not int or revision < 1:
            return f"The revision {revision!r} is not an integer of 1 or more."
        data = header.get("data") or []
        ids = [entry.get("id") for entry in data if isinstance(entry, dict)]
        if len(set(ids)) != len(ids):
            return "Two data arrays have the same id."
        for entry in data:
            problem = _check_data(
                object_type, entry, arrays, buffers, header.get("params"), header.get("log_sets")
            )
            if problem is not None:
                return problem
        for group in header.get("property_groups") or []:
            problem = _check_group(group, ids)
            if problem is not None:
                return problem
        for name, array in buffers.items():
            if array.size * array.itemsize != array.nbytes:
                return f"The buffer {name!r} does not match its declared size."
        size = sum(int(array.nbytes) for array in buffers.values())
        if size > MAX_OBJECT_BYTES:
            return f"The object is {size} bytes. The limit is {MAX_OBJECT_BYTES} bytes."
        return None

    def tool_frame(self, data: bytes) -> list[Frame]:
        """Apply one tool-client frame and return the reply frames.

        Parameters
        ----------
        data : bytes
            The frame from the client.

        Returns
        -------
        list of Frame
            The frames that the viewer sends back.
        """
        try:
            frame = protocol.decode_frame(data)
        except ProtocolError as error:
            return [_error("invalid_frame", str(error))]
        with self._state_lock:
            self.tool_requests.append(frame)
        return self._tool_reply(frame)

    def _tool_reply(self, frame: Frame) -> list[Frame]:
        """Return the reply frames for one decoded tool-client frame.

        For ``object.fetch`` the list holds all reply frames. A real connection
        uses ``_stream_fetch`` instead, which can stop on ``command.cancel``.
        """
        kind = frame.header.get("type")
        request = frame.header.get("request")
        if self.protocol_minor >= 11 and kind == "commands.list":
            return [self._commands_list(request)]
        if self.protocol_minor >= 11 and kind == "command.run":
            return [self._command_run(frame)]
        if self.protocol_minor >= 18 and kind == "object.fetch":
            return list(self._fetch_frames(frame))
        if self.protocol_minor >= 18 and kind in ("selection.set", "selection.get"):
            if self.kind == "headless":
                return [_request_error(request, "needs_window", "The viewer has no window.")]
            if kind == "selection.get":
                return [self._selection_get(request)]
            return [self._selection_set(frame)]
        return [_error("unknown_message", f"The message type {kind!r} is not known.")]

    def _commands_list(self, request: object) -> Frame:
        """Return the ``commands`` frame for ``commands.list``."""
        entries = []
        for name, (description, needs_window, fields) in _COMMANDS.items():
            properties = {key: dict(_INPUT_SCHEMAS[kind]) for key, (kind, _) in fields.items()}
            if name == "colour.set_map":
                properties["name"] = (
                    {"type": "string", "enum": list(_COLOUR_MAPS)}
                    if self.protocol_minor >= 18
                    else {"type": "string"}
                )
            entry: dict[str, Any] = {
                "name": name,
                "description": description,
                "needs_window": needs_window,
                "tool_available": True,
                "unavailable_reason": None,
                "input_schema": {
                    "type": "object",
                    "properties": properties,
                    "required": [key for key, (_, required) in fields.items() if required],
                    "additionalProperties": False,
                },
                "result_schema": {
                    "type": "object",
                    "properties": {
                        "status": {"enum": ["done", "refused"]},
                        "text": {"type": "string"},
                    },
                    "required": ["status"],
                },
            }
            if name == "filter.edit" and self.protocol_minor >= 20:
                entry["long"] = True
            entries.append(entry)
        header = {"type": "commands", "request": request, "viewer": self.kind, "commands": entries}
        return Frame(header=header, buffers={})

    def _wait_if_busy(self, command: str) -> bool:
        """Wait for the busy time to end. Return ``False`` if the command must get ``busy``.

        ``colour.set_scale`` does not wait. It gets ``busy`` at once. Other commands wait
        at most ``busy_wait`` seconds.
        """
        with self._state_lock:
            remaining = self._busy_until - time.monotonic()
        if remaining <= 0:
            return True
        if command == "colour.set_scale":
            return False
        time.sleep(min(remaining, self.busy_wait))
        return remaining <= self.busy_wait

    def _command_problem(self, frame: Frame) -> Frame | None:
        """Check the input and the object of a command. Return a ``request.error`` or ``None``."""
        request = frame.header.get("request")
        command = frame.header.get("command")
        spec = _COMMANDS.get(command) if isinstance(command, str) else None
        if spec is None:
            return _request_error(
                request, "invalid_input", f"The command {command!r} is not known."
            )
        _, needs_window, fields = spec
        given = frame.header.get("input")
        bad = _check_input(fields, given)
        if bad is None and command == "filter.edit":
            bad = _check_edit(given["edit"])
        if bad is not None:
            return _request_error(
                request, "invalid_input", f"The input field {bad!r} is not valid.", bad
            )
        id_field = "object" if "object" in fields else "id"
        with self._state_lock:
            known = given[id_field] in self.objects
        if not known:
            return _request_error(
                request,
                "invalid_input",
                f"The object {given[id_field]!r} is not in the scene.",
                id_field,
            )
        if needs_window and self.kind == "headless":
            return _request_error(request, "needs_window", "The viewer has no window.")
        return None

    def _command_run(self, frame: Frame) -> Frame:
        """Check and run one command. Return the ``command.result`` or ``request.error`` frame."""
        problem = self._command_problem(frame)
        if problem is not None:
            return problem
        request = frame.header.get("request")
        command = frame.header.get("command")
        given = frame.header["input"]
        cancel = command == "filter.edit" and given["edit"]["kind"] == "cancel"
        if not cancel and not self._wait_if_busy(command):
            return _request_error(request, "busy", "The scene is busy with a filter change.")
        if command == "filter.edit":
            result = self._filter_edit(given["object"], given["edit"])
        else:
            result = self._colour_command(command, given)
        header = {"type": "command.result", "request": request, "result": result}
        return Frame(header=header, buffers={})

    def _filter_edit(self, object_id: str, edit: dict[str, Any]) -> dict[str, Any]:
        """Apply one filter edit and return the result dict.

        ``add`` makes an off filter. Its kind is ``"values"`` for a ``referenced`` or
        ``boolean`` array, and ``"range"`` for other arrays. The ``set_*`` edits change
        a filter that exists. ``remove`` deletes a filter. ``cancel`` ends the busy time.
        An edit that changes a filter starts the busy time.
        """
        kind = edit["kind"]
        if kind == "cancel":
            with self._state_lock:
                self._busy_until = 0.0
            return {"status": "done"}
        array = edit["array"]
        with self._state_lock:
            entries = {
                e.get("name"): e
                for e in self.objects[object_id].header.get("data") or []
                if isinstance(e, dict)
            }
            entry = entries.get(array)
            filters = self._filters.setdefault(object_id, {})
            if kind == "add":
                if entry is None:
                    return {"status": "refused", "text": f"The object has no array {array!r}."}
                if entry.get("data_type") == "text" or entry.get("association") == "object":
                    return {
                        "status": "refused",
                        "text": f"The array {array!r} cannot have a filter.",
                    }
                by_values = entry.get("data_type") in _FIXED_COLOUR_TYPES
                filters[array] = {
                    "kind": "values" if by_values else "range",
                    "on": False,
                    "minimum": None,
                    "maximum": None,
                    "keep_nan": False,
                    "values": [],
                    "exclude": False,
                }
            elif kind == "remove":
                if filters.pop(array, None) is None:
                    return {"status": "done"}
            else:
                state = filters.get(array)
                if state is None:
                    return {"status": "refused", "text": f"The array {array!r} has no filter."}
                if kind == "set_range":
                    if state["kind"] != "range":
                        return {
                            "status": "refused",
                            "text": f"The array {array!r} has no range filter.",
                        }
                    try:
                        state["minimum"] = float(edit["min"])
                        state["maximum"] = float(edit["max"])
                    except ValueError:
                        return {"status": "refused", "text": "The range is not a number."}
                    state["keep_nan"] = edit["keep_nan"]
                elif kind == "set_values":
                    if state["kind"] != "values":
                        return {
                            "status": "refused",
                            "text": f"The array {array!r} has no values filter.",
                        }
                    state["values"] = list(edit["values"])
                elif kind == "set_on":
                    state["on"] = edit["on"]
                else:
                    state["exclude"] = edit["exclude"]
            if self.protocol_minor < 20:
                self._busy_until = time.monotonic() + self.filter_delay
        return {"status": "done"}

    def _filter_counts(self, object_id: str, array: str) -> tuple[int, int]:
        """Return ``(passing, total)`` for the array ``array`` of an object.

        The count uses the filters that are on and that have the association of the
        array. A range filter keeps ``minimum <= v <= maximum``, and NaN only if
        ``keep_nan`` is set. A values filter keeps the listed codes, or all other
        codes if ``exclude`` is set.
        """
        with self._state_lock:
            stored = self.objects[object_id]
            entries = {
                e.get("name"): e for e in stored.header.get("data") or [] if isinstance(e, dict)
            }
            filters = copy.deepcopy(self._filters.get(object_id, {}))
            buffers = stored.buffers
        entry = entries[array]
        total = int(buffers[entry["buffer"]].shape[0])
        keep = np.ones(total, dtype=bool)
        for name, state in filters.items():
            other = entries.get(name)
            if not state["on"] or other is None:
                continue
            if other.get("association") != entry.get("association"):
                continue
            values = buffers[other["buffer"]]
            if state["kind"] == "range":
                passed = (values >= state["minimum"]) & (values <= state["maximum"])
                if state["keep_nan"]:
                    passed = passed | np.isnan(values)
            else:
                passed = np.isin(values, state["values"])
                if state["exclude"]:
                    passed = ~passed
            keep &= passed
        return int(keep.sum()), total

    def _build(self, conn: socket.socket, reader: _Reader, request: object) -> str:
        """Send ``command.progress`` frames for one filter build.

        The method sends a frame about every 0.1 s and reads a possible
        ``command.cancel`` between them. It returns ``"cancelled"`` on cancel,
        ``"stopped"`` when ``filter_stop_after`` ends the wait, and ``"done"`` at
        the end of the build.
        """
        version = protocol.TOOL_PROTOCOL_VERSION
        stop = self.filter_stop_after
        if stop is not None and stop >= self.filter_delay:
            stop = None
        started = time.monotonic()
        fraction = 0.0
        while True:
            progress = {
                "type": "command.progress",
                "request": request,
                "fraction": fraction,
                "text": "filter.edit",
            }
            _send(conn, _encode(Frame(header=progress, buffers={}), version))
            if fraction >= 1.0:
                return "done"
            end = self.filter_delay if stop is None else stop
            wait = min(0.1, end - (time.monotonic() - started))
            if wait > 0:
                time.sleep(wait)
            if self._cancelled(conn, reader, request):
                return "cancelled"
            elapsed = time.monotonic() - started
            if stop is not None and elapsed >= stop:
                return "stopped"
            fraction = min(1.0, elapsed / self.filter_delay)

    def _filter_request(self, data: bytes) -> Frame | None:
        """Return the frame if ``data`` is a ``filter.edit`` that this viewer builds in steps.

        The method also keeps the frame in ``tool_requests``. It returns ``None`` below
        protocol 1.20 and for other messages.
        """
        if self.protocol_minor < 20:
            return None
        try:
            frame = protocol.decode_frame(data)
        except ProtocolError:
            return None
        if frame.header.get("type") != "command.run":
            return None
        if frame.header.get("command") != "filter.edit":
            return None
        with self._state_lock:
            self.tool_requests.append(frame)
        return frame

    def _stream_filter(self, conn: socket.socket, reader: _Reader, frame: Frame) -> None:
        """Run one ``filter.edit`` and send ``command.progress`` frames during a build.

        A build runs when the edit is ``set_range``, ``set_values``, ``set_exclude`` or
        ``set_on``, the filter is on after the edit and ``filter_delay`` is above 0. On
        ``command.cancel``, the method restores the filters of the object and sends
        ``cancelled``. If ``filter_stop_after`` is shorter than the build, the result is
        ``{"status": "done", "state": "building"}`` at that time. At the end of a build,
        and for an edit with no build, the result has ``passing`` and ``total``.
        """
        version = protocol.TOOL_PROTOCOL_VERSION
        request = frame.header.get("request")
        problem = self._command_problem(frame)
        if problem is not None:
            _send(conn, _encode(problem, version))
            return
        given = frame.header["input"]
        object_id = given["object"]
        edit = given["edit"]
        with self._state_lock:
            before = copy.deepcopy(self._filters.get(object_id))
        result = self._filter_edit(object_id, edit)
        on = False
        if result["status"] == "done" and edit["kind"].startswith("set_"):
            with self._state_lock:
                on = self._filters[object_id][edit["array"]]["on"]
        # ``set_on`` always gives the count. Other edits give it only if the filter is on.
        if result["status"] == "done" and (on or edit["kind"] == "set_on"):
            array = edit["array"]
            outcome = "done"
            if on and self.filter_delay > 0:
                outcome = self._build(conn, reader, request)
            if outcome == "cancelled":
                with self._state_lock:
                    if before is None:
                        self._filters.pop(object_id, None)
                    else:
                        self._filters[object_id] = before
                error = _request_error(request, "cancelled", "The client cancelled the request.")
                _send(conn, _encode(error, version))
                return
            if outcome == "stopped":
                result = {"status": "done", "state": "building"}
            else:
                passing, total = self._filter_counts(object_id, array)
                result = {"status": "done", "passing": passing, "total": total}
        header = {"type": "command.result", "request": request, "result": result}
        _send(conn, _encode(Frame(header=header, buffers={}), version))

    def _colour_command(self, command: str, given: dict[str, Any]) -> dict[str, Any]:
        """Run one colour command and return its result dict."""
        object_id = given["object"] if command == "colour.set_scale" else given["id"]
        with self._state_lock:
            stored = self.objects[object_id]
            state = self._colouring.setdefault(
                object_id,
                {
                    "array": None,
                    "colour_map": None,
                    "range": None,
                    "scale": "linear",
                    "plain": None,
                },
            )
            entries = {
                e.get("name"): e for e in stored.header.get("data") or [] if isinstance(e, dict)
            }
            if command == "colour.set_active":
                array = given.get("array")
                if array is not None and array not in entries:
                    return {
                        "status": "refused",
                        "text": f"The object {object_id!r} has no array {array!r}.",
                    }
                state["array"] = array
            elif command == "colour.clear_range":
                state["range"] = None
            elif command == "object.set_colour":
                colour = given["colour"]
                state["plain"] = None if colour is None else list(colour)
            else:
                active = entries.get(state["array"])
                if active is not None and active.get("data_type") in _FIXED_COLOUR_TYPES:
                    return {
                        "status": "refused",
                        "text": f"The array {state['array']!r} has a fixed colouring.",
                    }
                if command == "colour.set_map":
                    if self.protocol_minor >= 18 and given["name"] not in _COLOUR_MAPS:
                        return {
                            "status": "refused",
                            "text": f"The colour map {given['name']!r} is not known.",
                        }
                    state["colour_map"] = given["name"]
                elif command == "colour.set_range":
                    try:
                        state["range"] = (float(given["min"]), float(given["max"]))
                    except ValueError:
                        return {"status": "refused", "text": "The range is not a number."}
                else:
                    state["scale"] = given["scale"]
        return {"status": "done"}

    def _selection_get(self, request: object) -> Frame:
        """Return the ``selection.result`` frame for ``selection.get``."""
        header: dict[str, Any] = {"type": "selection.result", "request": request, "entries": []}
        buffers: dict[str, np.ndarray] = {}
        with self._state_lock:
            current = self._selection
            if current is not None:
                header["object"] = current["object"]
                header["revision"] = current["revision"]
                for kind in ("cells", "nodes"):
                    values = current[kind]
                    if values is not None:
                        header["entries"].append(
                            {"kind": kind, "count": int(values.size), "indices": kind}
                        )
                        buffers[kind] = np.asarray(values, dtype=np.uint32)
        return Frame(header=header, buffers=buffers)

    def _selection_set(self, frame: Frame) -> Frame:
        """Check and keep a selection. Return the reply to ``selection.set``."""
        header = frame.header
        request = header.get("request")
        object_id = header.get("object")
        entries = header.get("entries")
        if not isinstance(entries, list):
            return _request_error(request, "invalid_input", "The entries are missing.", "entries")
        kinds: dict[str, np.ndarray] = {}
        for item in entries:
            name = item.get("indices") if isinstance(item, dict) else None
            array = frame.buffers.get(name) if isinstance(name, str) else None
            if array is None or array.dtype != np.dtype("uint32") or array.ndim != 1:
                return _request_error(
                    request, "invalid_input", "An entry needs a u32 buffer of indices.", "indices"
                )
            kinds[item.get("kind")] = array.astype(np.int64)
        with self._state_lock:
            problem = self._selection_problem(object_id, header.get("revision"), kinds)
            if problem is not None:
                return _request_error(request, "invalid_input", problem[1], problem[0])
            self._store_selection(object_id, kinds)
            revision = self.objects[object_id].header.get("revision", 1)
        result: dict[str, Any] = {
            "type": "selection.result",
            "request": request,
            "object": object_id,
            "revision": revision,
            "entries": [
                {"kind": kind, "count": int(np.unique(values).size)}
                for kind, values in kinds.items()
            ],
        }
        return Frame(header=result, buffers={})

    def _fetch_frames(self, frame: Frame) -> Iterator[Frame]:
        """Yield the reply frames of ``object.fetch``.

        The frames are the ``object.header`` frame, the ``object.part`` frames of each
        array in order (each followed by ``command.progress`` if there is more than one
        part), and the ``object.end`` frame. An unknown id gives one ``request.error``.
        """
        request = frame.header.get("request")
        object_id = frame.header.get("id")
        with self._state_lock:
            stored = self.objects.get(object_id) if isinstance(object_id, str) else None
            if stored is None:
                header = None
                buffers: dict[str, np.ndarray] = {}
            else:
                header = dict(stored.header, buffers=[])
                header["data"] = [
                    {"offsets": None, "value_map": None, "colour_map": None, **entry}
                    for entry in stored.header.get("data") or []
                ]
                buffers = dict(stored.buffers)
        if header is None:
            yield _request_error(
                request, "invalid_input", f"The object {object_id!r} is not in the scene.", "id"
            )
            return
        step = max(8, self.part_bytes // 8 * 8)
        raw = {
            name: np.ascontiguousarray(array)
            .astype(array.dtype.newbyteorder("<"), copy=False)
            .tobytes()
            for name, array in buffers.items()
        }
        total = sum(-(-len(data) // step) for data in raw.values())
        arrays = [
            {
                "name": name,
                "dtype": _PROTOCOL_DTYPES[array.dtype.name],
                "shape": [int(n) for n in array.shape],
            }
            for name, array in buffers.items()
        ]
        yield Frame(
            header={
                "type": "object.header",
                "request": request,
                "origin": list(self.origin),
                "crs": self.crs,
                "object": header,
                "arrays": arrays,
            },
            buffers={},
        )
        sent = 0
        for name, data in raw.items():
            for offset in range(0, len(data), step):
                chunk = data[offset : offset + step]
                sent += 1
                yield Frame(
                    header={
                        "type": "object.part",
                        "request": request,
                        "array": name,
                        "offset": offset,
                        "bytes": len(chunk),
                    },
                    buffers={"part": np.frombuffer(chunk, dtype=np.uint8)},
                )
                if total > 1:
                    progress = {
                        "type": "command.progress",
                        "request": request,
                        "fraction": sent / total,
                        "text": f"Part {sent} of {total}",
                    }
                    yield Frame(header=progress, buffers={})
        yield Frame(header={"type": "object.end", "request": request, "parts": total}, buffers={})

    def _stream_fetch(self, conn: socket.socket, reader: _Reader, frame: Frame) -> None:
        """Send the reply to ``object.fetch`` and read between the parts.

        After the header and after each part, the method waits at most 0.001 s for a
        message (an object without parts has the header only). If it is a
        ``command.cancel`` with the same request, it sends ``cancelled`` and stops.
        """
        request = frame.header.get("request")
        for reply in self._fetch_frames(frame):
            _send(conn, _encode(reply, protocol.TOOL_PROTOCOL_VERSION))
            if reply.header["type"] not in ("object.header", "object.part"):
                continue
            if self.fetch_delay:
                time.sleep(self.fetch_delay)
            if self._cancelled(conn, reader, request):
                error = _request_error(request, "cancelled", "The client cancelled the request.")
                _send(conn, _encode(error, protocol.TOOL_PROTOCOL_VERSION))
                return

    @staticmethod
    def _cancelled(conn: socket.socket, reader: _Reader, request: object) -> bool:
        """Return whether a ``command.cancel`` for ``request`` arrives within 0.001 s.

        A message that is not a ``command.cancel`` is kept in ``reader.deferred``.
        The viewer handles it after the current work, as a real viewer does. A
        ``command.cancel`` for another request is dropped.
        """
        previous = conn.gettimeout()
        conn.settimeout(0.001)
        try:
            data = reader.message()
        except (TimeoutError, BlockingIOError):
            return False
        finally:
            conn.settimeout(previous)
        if data is None:
            raise _Closed
        try:
            header = protocol.decode_frame(data).header
        except ProtocolError:
            reader.deferred.append(data)
            return False
        if header.get("type") != "command.cancel":
            reader.deferred.append(data)
            return False
        return header.get("request") == request

    # -- server -------------------------------------------------------------

    def _serve(self) -> None:
        """Accept connections, each in its own thread."""
        while True:
            try:
                conn, _ = self._server.accept()
            except OSError:
                return
            with self._lock:
                self._connections.add(conn)
            threading.Thread(target=self._connection, args=(conn,), daemon=True).start()

    def _connection(self, conn: socket.socket) -> None:
        """Run one connection until it ends."""
        is_data = False
        try:
            opened = self._handshake(conn)
            if opened is None:
                return
            path, reader = opened
            if path.split("?")[0] == "/tool":
                self._run(conn, reader, tool=True)
                return
            with self._lock:
                busy = self._data_connections > 0
                if not busy:
                    self._data_connections += 1
                    is_data = True
            if busy:
                refusal = _error("busy", "Another client is connected.")
                _send(conn, _encode(refusal))
                conn.sendall(bytes([0x88, 0]))
                self._drain(conn)
                return
            conn.settimeout(self.heartbeat_timeout)
            self._run(conn, reader, tool=False)
        except (_Closed, OSError):
            pass
        finally:
            if is_data:
                with self._lock:
                    self._data_connections -= 1
            with self._lock:
                self._connections.discard(conn)
            try:
                conn.close()
            except OSError:
                pass

    @staticmethod
    def _drain(conn: socket.socket) -> None:
        """Stop sending, then read until the peer closes or 0.5 s pass.

        This stops a reset from the system that would destroy the last frames.
        """
        try:
            conn.shutdown(socket.SHUT_WR)
            conn.settimeout(0.5)
            while conn.recv(65536):
                pass
        except OSError:
            pass

    def _handshake(self, conn: socket.socket) -> tuple[str, _Reader] | None:
        """Read the HTTP upgrade request. Return the path and a reader, or ``None``."""
        data = b""
        while b"\r\n\r\n" not in data:
            chunk = conn.recv(4096)
            if not chunk:
                return None
            data += chunk
        head, _, rest = data.partition(b"\r\n\r\n")
        lines = head.decode("latin-1").split("\r\n")
        parts = lines[0].split(" ")
        path = parts[1] if len(parts) > 1 else "/"
        headers = {
            k.strip().lower(): v.strip() for k, _, v in (line.partition(":") for line in lines[1:])
        }
        if headers.get("authorization") != f"Bearer {self.token}":
            conn.sendall(b"HTTP/1.1 401 Unauthorized\r\nContent-Length: 0\r\n\r\n")
            return None
        key = headers.get("sec-websocket-key")
        if key is None:
            conn.sendall(b"HTTP/1.1 400 Bad Request\r\nContent-Length: 0\r\n\r\n")
            return None
        accept = base64.b64encode(hashlib.sha1(key.encode() + _GUID).digest()).decode()
        conn.sendall(
            (
                "HTTP/1.1 101 Switching Protocols\r\nUpgrade: websocket\r\n"
                f"Connection: Upgrade\r\nSec-WebSocket-Accept: {accept}\r\n\r\n"
            ).encode()
        )
        return path, _Reader(conn, rest)

    def _fetch_request(self, data: bytes) -> Frame | None:
        """Return the frame if ``data`` is an ``object.fetch`` that this viewer supports.

        The method also keeps the frame in ``tool_requests``.
        """
        if self.protocol_minor < 18:
            return None
        try:
            frame = protocol.decode_frame(data)
        except ProtocolError:
            return None
        if frame.header.get("type") != "object.fetch":
            return None
        with self._state_lock:
            self.tool_requests.append(frame)
        return frame

    def _run(self, conn: socket.socket, reader: _Reader, tool: bool) -> None:
        """Read messages and send the replies.

        The method handles the messages in ``reader.deferred`` first, in order.
        """
        while True:
            data = reader.deferred.pop(0) if reader.deferred else reader.message()
            if data is None:
                return
            if tool:
                fetch = self._fetch_request(data)
                if fetch is not None:
                    self._stream_fetch(conn, reader, fetch)
                    continue
                edit = self._filter_request(data)
                if edit is not None:
                    self._stream_filter(conn, reader, edit)
                    continue
                for reply in self.tool_frame(data):
                    _send(conn, _encode(reply, protocol.TOOL_PROTOCOL_VERSION))
            else:
                for reply in self.accept_frame(data):
                    _send(conn, _encode(reply))
            if self.read_delay:
                time.sleep(self.read_delay)
