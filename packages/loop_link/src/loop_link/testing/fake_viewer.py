"""A fake viewer for the tests: a WebSocket server on 127.0.0.1.

The fake viewer does what the link of the real viewer does. It checks the bearer
token, answers ``hello`` with ``hello`` and ``scene.inventory``, and keeps the
objects that it receives. It has its own small server-side WebSocket code.
"""

from __future__ import annotations

import base64
import hashlib
import json
import math
import os
import socket
import struct
import threading
import time
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


def _encode(frame: Frame) -> bytes:
    """Return the bytes of a frame."""
    return protocol.encode_frame(frame.header, list(frame.buffers.items()))


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
    """Read masked client frames from a socket."""

    def __init__(self, conn: socket.socket, buffer: bytes = b"") -> None:
        self._conn = conn
        self._buffer = buffer

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

    Attributes
    ----------
    port : int
        The port of the server.
    objects : dict of str to StoredObject
        The objects in the scene.
    received : list of Frame
        All frames that the viewer received on all connections.
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
    ) -> None:
        self.token = token
        self.crs = crs
        self.origin = origin
        self.kind = kind
        self.heartbeat_timeout = heartbeat_timeout
        self.reply_errors = dict(reply_errors or {})
        self.read_delay = read_delay
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

    def __enter__(self) -> FakeViewer:  # noqa: PYI034
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
        with self._state_lock:
            self.objects[object_id] = StoredObject(header=header, buffers=dict(frame.buffers))
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

    def _run(self, conn: socket.socket, reader: _Reader, tool: bool) -> None:
        """Read messages and send the replies."""
        while True:
            data = reader.message()
            if data is None:
                return
            if tool:
                replies = [_error("unknown_message", "The tool endpoint has no messages.")]
            else:
                replies = self.accept_frame(data)
            for reply in replies:
                _send(conn, _encode(reply))
            if self.read_delay:
                time.sleep(self.read_delay)
