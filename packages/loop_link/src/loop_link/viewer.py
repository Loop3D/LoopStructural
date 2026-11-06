"""The ``Viewer`` class: one connection from Python to the viewer."""

from __future__ import annotations

import threading
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import Any, get_args

from loop_link import discovery, encode, protocol
from loop_link._websocket import WebSocketClient, WebSocketClosed, WebSocketError
from loop_link.errors import (
    BusyError,
    InvalidInputError,
    MissingFeatureError,
    Refusal,
    ViewerClosedError,
    ViewerRefusedError,
)
from loop_link.objects import SceneObject

_ACCEPTED_INPUTS = (
    "Accepted inputs: loop_link objects (Points, Vectors, Lines, Surface, BlockModel, "
    "TetraMesh, Drillhole); pyvista meshes and grids; LoopStructural Surface and "
    "StructuredGrid; pandas tables for drillholes "
    "(use loop_link.adapters.pandas.drillholes_from_pandas)."
)
_LINK_ERRORS = (OSError, WebSocketClosed, WebSocketError)


@dataclass(frozen=True)
class InventoryEntry:
    """One object in the scene of the viewer.

    Parameters
    ----------
    id : str
        The identifier of the object.
    object_type : str
        The type of the object, as the viewer names it.
    revision : int
        The revision of the object in the viewer.
    """

    id: str
    object_type: str
    revision: int


class Viewer:
    """A connection to the viewer.

    Parameters
    ----------
    port : int, optional
        The port of the viewer. The default is to search for it.
    token : str, optional
        The access token. The default is to read it from the environment or the settings.
    host : str
        The host of the viewer.
    crs : str
        The coordinate reference system that this client proposes.
    origin : sequence of float
        The scene origin that this client proposes, three values.
    settings_dir : pathlib.Path, optional
        The settings folder of the viewer.
    timeout : float
        The time in seconds to wait for the connection and the first replies.
    heartbeat_interval : float
        The time in seconds between keepalive checks.
    """

    def __init__(
        self,
        *,
        port: int | None = None,
        token: str | None = None,
        host: str = "127.0.0.1",
        crs: str = "",
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        settings_dir: Path | None = None,
        timeout: float = 5.0,
        heartbeat_interval: float = 5.0,
    ) -> None:
        self._port = port
        self._token = token
        self._host = host
        self._crs = crs
        self._proposed_origin = _triple(origin)
        self._settings_dir = settings_dir
        self._timeout = timeout
        self._heartbeat_interval = heartbeat_interval

        self._lock = threading.Lock()
        self._client: WebSocketClient | None = None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._dead = False
        self._last_send = 0.0
        self._seq = 0

        self._revisions: dict[str, int] = {}
        self._inventory: list[InventoryEntry] = []
        self._inventory_revisions: dict[str, int] = {}
        self._refusals: list[Refusal] = []
        self._scene_origin = self._proposed_origin
        self._scene_crs = crs

    # -- properties -----------------------------------------------------------

    @property
    def revisions(self) -> Mapping[str, int]:
        """The last revision that this client used for each id."""
        return dict(self._revisions)

    @property
    def scene_origin(self) -> tuple[float, float, float]:
        """The scene origin that all frames use."""
        return self._scene_origin

    @property
    def scene_crs(self) -> str:
        """The coordinate reference system of the scene, as the viewer reported it."""
        return self._scene_crs

    @property
    def connected(self) -> bool:
        """True while the connection is open."""
        return self._client is not None and not self._dead

    # -- context manager ------------------------------------------------------

    def __enter__(self) -> Viewer:  # noqa: PYI034
        """Connect and return the viewer."""
        self.connect()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Close the connection."""
        self.close()

    # -- connection -----------------------------------------------------------

    def connect(self) -> None:
        """Connect to the viewer and read its ``hello`` and inventory.

        Raises
        ------
        loop_link.errors.NoViewerError
            No viewer was found.
        loop_link.errors.TokenError
            The viewer refused the token.
        loop_link.errors.BusyError
            Another client owns the scene.
        loop_link.errors.MissingFeatureError
            The viewer does not support this protocol version.
        loop_link.errors.ViewerRefusedError
            The viewer sent another error.
        loop_link.errors.ViewerClosedError
            The viewer closed the connection before it finished the hello.
        """
        from loop_link import __version__

        if self.connected:
            return
        self._stop_thread()
        self._drop_client()
        address = discovery.find_viewer(
            port=self._port,
            token=self._token,
            host=self._host,
            settings=self._settings_dir,
        )
        client = WebSocketClient.connect(
            address.host, address.port, token=address.token, timeout=self._timeout, path="/"
        )
        try:
            client.send_binary(encode.hello_frame(self._crs, self._proposed_origin, __version__))
            hello, entries = self._read_hello(client)
        except _LINK_ERRORS as error:
            _close_client(client)
            raise ViewerClosedError from error
        except BaseException:
            _close_client(client)
            raise
        origin = hello.get("origin")
        self._scene_origin = _triple(origin) if origin is not None else self._proposed_origin
        self._scene_crs = str(hello.get("crs") or "")
        self._inventory = entries
        self._inventory_revisions = {entry.id: entry.revision for entry in entries}
        self._refusals = []
        self._client = client
        self._dead = False
        self._last_send = time.monotonic()
        self._stop = threading.Event()
        self._thread = threading.Thread(
            target=self._heartbeat_loop,
            args=(self._stop,),
            name="loop-link-heartbeat",
            daemon=True,
        )
        self._thread.start()

    def _read_hello(self, client: WebSocketClient) -> tuple[dict[str, Any], list[InventoryEntry]]:
        """Read frames until the ``hello`` and the inventory arrive."""
        hello: dict[str, Any] | None = None
        entries: list[InventoryEntry] | None = None
        deadline = time.monotonic() + self._timeout
        while hello is None or entries is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise ViewerClosedError(
                    "The viewer did not answer the hello in time. Check the viewer and retry."
                )
            data = client.receive(min(remaining, 0.5))
            if data is None:
                continue
            header = protocol.decode_frame(data).header
            kind = header.get("type")
            if kind == "error":
                _raise_hello_error(header)
            elif kind == "hello":
                hello = header
            elif kind == "scene.inventory":
                entries = [
                    InventoryEntry(
                        id=str(item["id"]),
                        object_type=str(item.get("object_type", "")),
                        revision=int(item.get("revision", 1)),
                    )
                    for item in header.get("entries", [])
                ]
        return hello, entries

    # -- sending --------------------------------------------------------------

    def send(
        self,
        obj: SceneObject,
        *,
        id: str,
        name: str | None = None,
        parent: str | None = None,
    ) -> int:
        """Send one object to the viewer.

        Parameters
        ----------
        obj : SceneObject
            The object to send.
        id : str
            The identifier of the object. The same id replaces the object.
        name : str, optional
            The display name. The default is the id.
        parent : str, optional
            The identifier of a group that holds the object.

        Returns
        -------
        int
            The revision that the frame used.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The input type is not accepted.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        if not isinstance(obj, get_args(SceneObject)):
            raise InvalidInputError(_ACCEPTED_INPUTS)
        self._ensure_connected()
        revision = self._next_revision(id)
        data = encode.encode_object(
            obj,
            id=id,
            name=id if name is None else name,
            parent=parent,
            revision=revision,
            origin=self._scene_origin,
        )
        self._send_frame(data)
        return revision

    def send_group(self, id: str, name: str | None = None, parent: str | None = None) -> int:
        """Send a group.

        Parameters
        ----------
        id : str
            The identifier of the group.
        name : str, optional
            The display name. The default is the id.
        parent : str, optional
            The identifier of a parent group.

        Returns
        -------
        int
            The revision that the frame used.
        """
        self._ensure_connected()
        revision = self._next_revision(id)
        self._send_frame(encode.group_frame(id, id if name is None else name, parent, revision))
        return revision

    def remove(self, id: str) -> None:
        """Remove one object from the scene. The revision count of the id stays.

        Parameters
        ----------
        id : str
            The identifier of the object.
        """
        self._ensure_connected()
        self._send_frame(encode.remove_frame(id))

    def clear(self) -> None:
        """Remove all objects from the scene."""
        self._ensure_connected()
        self._send_frame(encode.clear_frame())

    def inventory(self) -> list[InventoryEntry]:
        """Return the inventory that the viewer sent at hello time.

        Returns
        -------
        list of InventoryEntry
            The entries. The list is not updated after the hello.
        """
        return list(self._inventory)

    def flush(self, settle: float = 0.5, *, raise_on_refusal: bool = True) -> list[Refusal]:
        """Wait for late replies and return the refusals.

        Parameters
        ----------
        settle : float
            The time in seconds to wait for late frames.
        raise_on_refusal : bool
            If true, raise when the viewer refused an object.

        Returns
        -------
        list of Refusal
            The refusals since the last flush.

        Raises
        ------
        loop_link.errors.ViewerRefusedError
            The viewer refused one or more objects and ``raise_on_refusal`` is true.
        """
        self._ensure_connected()
        deadline = time.monotonic() + settle
        while True:
            self._drain(max(0.001, min(0.05, deadline - time.monotonic())))
            if time.monotonic() >= deadline:
                break
        refusals, self._refusals = self._refusals, []
        if refusals and raise_on_refusal:
            raise ViewerRefusedError(refusals)
        return refusals

    # -- closing --------------------------------------------------------------

    def release_data(self) -> None:
        """Close the connection so another client can own the scene."""
        self.close()

    def close(self) -> None:
        """Close the connection. The revision numbers stay."""
        self._stop_thread()
        self._drop_client()
        self._dead = True

    # -- helpers --------------------------------------------------------------

    def _ensure_connected(self) -> None:
        if self._dead:
            raise ViewerClosedError
        if self._client is None:
            self.connect()

    def _drop_client(self) -> None:
        """Close the socket, if there is one."""
        with self._lock:
            client, self._client = self._client, None
        if client is not None:
            _close_client(client)

    def _next_revision(self, id: str) -> int:
        revision = max(self._revisions.get(id, 0), self._inventory_revisions.get(id, 0)) + 1
        self._revisions[id] = revision
        return revision

    def _send_frame(self, data: bytes) -> None:
        with self._lock:
            client = self._client
            if client is None:
                raise ViewerClosedError
            try:
                client.send_binary(data)
                self._last_send = time.monotonic()
                self._drain_locked(client, 0.001)
            except _LINK_ERRORS as error:
                self._dead = True
                raise ViewerClosedError from error

    def _drain(self, wait: float) -> None:
        with self._lock:
            client = self._client
            if client is None:
                raise ViewerClosedError
            try:
                self._drain_locked(client, wait)
            except _LINK_ERRORS as error:
                self._dead = True
                raise ViewerClosedError from error

    def _drain_locked(self, client: WebSocketClient, wait: float) -> None:
        """Read all waiting frames. The caller holds the lock."""
        while True:
            data = client.receive(wait)
            if data is None:
                return
            header = protocol.decode_frame(data).header
            if header.get("type") == "error" and header.get("object") is not None:
                self._refusals.append(
                    Refusal(
                        str(header["object"]),
                        str(header.get("code", "")),
                        str(header.get("message", "")),
                    )
                )
            wait = 0.001

    def _heartbeat_loop(self, stop: threading.Event) -> None:
        """Send a heartbeat when the connection is idle. Stop on a socket error."""
        interval = self._heartbeat_interval
        while not stop.wait(interval):
            if not self._lock.acquire(blocking=False):
                continue
            try:
                client = self._client
                if client is None:
                    return
                if time.monotonic() - self._last_send < interval:
                    continue
                self._seq += 1
                try:
                    client.send_binary(encode.heartbeat_frame(self._seq))
                    self._last_send = time.monotonic()
                except _LINK_ERRORS:
                    self._dead = True
                    return
            finally:
                self._lock.release()

    def _stop_thread(self) -> None:
        thread, self._thread = self._thread, None
        self._stop.set()
        if thread is not None and thread is not threading.current_thread():
            thread.join(timeout=2 * self._heartbeat_interval)


def _triple(values: Sequence[float]) -> tuple[float, float, float]:
    """Return three floats from a sequence of three numbers."""
    x, y, z = (float(v) for v in values)
    return (x, y, z)


def _close_client(client: WebSocketClient) -> None:
    """Close a client. An error on close does not matter."""
    try:
        client.close()
    except Exception:  # noqa: BLE001, S110
        pass


def _raise_hello_error(header: dict[str, Any]) -> None:
    """Raise the right error for an ``error`` frame during the hello."""
    code = header.get("code")
    message = str(header.get("message", ""))
    if code == "busy":
        raise BusyError(message or "The viewer is busy. Wait a moment and try again.")
    if code == "unsupported_version":
        raise MissingFeatureError("protocol version")
    raise ViewerRefusedError([Refusal(str(header.get("object", "")), str(code), message)])
