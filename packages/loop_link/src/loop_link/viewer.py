"""The ``Viewer`` class: one connection from Python to the viewer."""

from __future__ import annotations

import threading
import time
import weakref
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import Any

import numpy as np

from loop_link import discovery, encode, protocol
from loop_link._websocket import WebSocketClient, WebSocketClosed, WebSocketError
from loop_link.adapters import convert
from loop_link.commands import (
    ColourSettings,
    FilterResult,
    check_codes,
    check_number,
    check_rgb,
    check_text,
    format_number,
    list_colour_map_names,
    run_command,
)
from loop_link.errors import (
    BusyError,
    InvalidInputError,
    MissingFeatureError,
    NoViewerError,
    Refusal,
    ViewerClosedError,
    ViewerRefusedError,
)
from loop_link.fetch import FetchResult, fetch_object
from loop_link.protocol import ProtocolError
from loop_link.selection import Selection, SelectionWatcher, parse_selection_result
from loop_link.tool import ToolConnection

_LINK_ERRORS = (OSError, WebSocketClosed, WebSocketError)
_SELECTION_FEATURE = "selection and fetch (geoviewer protocol 1.18)"
_KIND_NAMES = {"cells": "cell", "nodes": "vertex"}
_INDEX_MAX = 2**32 - 1
_COLOUR_SCALES = ("linear", "log")


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
        self._shared: dict[str, Any] = {"client": None, "stop": self._stop}
        weakref.finalize(self, _release_shared, self._shared)
        self._thread: threading.Thread | None = None
        # _dead: the data connection failed. A later send opens a new connection.
        # _closed: close() was called. Only connect() opens the viewer again.
        self._dead = False
        self._closed = False
        self._address_lock = threading.Lock()
        self._address: discovery.ViewerAddress | None = None
        self._last_send = 0.0
        self._seq = 0
        self._tool_lock = threading.Lock()
        self._tool_connection: ToolConnection | None = None
        self._watchers: list[SelectionWatcher] = []
        self._colour_map_names: tuple[str, ...] | None = None

        self._revisions: dict[str, int] = {}
        self._inventory: dict[str, InventoryEntry] = {}
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
        client = self._client
        return client is not None and not client.closed and not self._dead

    @property
    def tool_connected(self) -> bool:
        """True while the tool connection is open."""
        tool = self._tool_connection
        return tool is not None and tool.connected and not self._closed

    # -- context manager ------------------------------------------------------

    def __enter__(self) -> Viewer:
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

    def connect(self) -> Viewer:
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

        Returns
        -------
        Viewer
            The viewer itself.
        """
        from loop_link import __version__

        if self.connected:
            return self
        self._closed = False
        self._stop_thread()
        self._drop_client()
        address = self._resolve_address()
        try:
            client = WebSocketClient.connect(
                address.host, address.port, token=address.token, timeout=self._timeout, path="/"
            )
        except NoViewerError:
            # The viewer at the known address is gone. Search again.
            self._forget_address()
            address = self._resolve_address()
            client = WebSocketClient.connect(
                address.host, address.port, token=address.token, timeout=self._timeout, path="/"
            )
        try:
            _read_early(client)
            client.send_binary(encode.hello_frame(self._crs, self._proposed_origin, __version__))
            hello, entries = self._read_hello(client)
            origin, crs = _parse_hello(hello)
        except _LINK_ERRORS as error:
            try:
                _read_early(client, 0.3)
            except _LINK_ERRORS:
                pass
            except BaseException:
                _close_client(client)
                raise
            _close_client(client)
            raise ViewerClosedError from error
        except BaseException:
            _close_client(client)
            raise
        self._scene_origin = origin if origin is not None else self._proposed_origin
        self._scene_crs = crs
        self._inventory = {entry.id: entry for entry in entries}
        self._inventory_revisions = {entry.id: entry.revision for entry in entries}
        self._refusals = []
        self._client = client
        self._shared["client"] = client
        self._dead = False
        self._last_send = time.monotonic()
        self._stop = threading.Event()
        self._shared["stop"] = self._stop
        self._thread = threading.Thread(
            target=_heartbeat_loop,
            args=(weakref.ref(self), self._stop, self._heartbeat_interval),
            name="loop-link-heartbeat",
            daemon=True,
        )
        self._thread.start()
        return self

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
                entries = _parse_inventory(header)
        return hello, entries

    # -- sending --------------------------------------------------------------

    def send(
        self,
        obj: Any,  # noqa: ANN401
        *,
        id: str,
        name: str | None = None,
        parent: str | None = None,
        colour: ColourSettings | None = None,
        **hints: Any,  # noqa: ANN401
    ) -> int:
        """Send one object to the viewer.

        If ``colour`` is given, the method waits for a late refusal of the object
        and then applies the colour settings with ``apply_colour``.

        Parameters
        ----------
        obj : Any
            The object to send: a loop_link object, or an input that an adapter converts.
        id : str
            The identifier of the object. The same id replaces the object.
        name : str, optional
            The display name. The default is the id.
        parent : str, optional
            The identifier of a group that holds the object.
        colour : loop_link.commands.ColourSettings, optional
            The colour settings to apply after the object is sent.
        **hints : Any
            Options for the adapter that converts ``obj``.

        Returns
        -------
        int
            The revision that the frame used.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The input type is not accepted, ``colour`` is not a ``ColourSettings``
            object, or the viewer refused a colour command. The object stays in
            the viewer.
        loop_link.errors.ViewerRefusedError
            The viewer refused this object and ``colour`` is given. Refusals of
            other objects stay for the next ``flush``.
        loop_link.errors.MissingFeatureError
            ``colour`` is given and the viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        if colour is not None and not isinstance(colour, ColourSettings):
            raise InvalidInputError("colour must be a ColourSettings object.")
        obj = convert(obj, **hints)
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
        self._revisions[id] = revision
        object_type = str(protocol.decode_header(data).get("object_type", ""))
        self._inventory[id] = InventoryEntry(id, object_type, revision)
        if colour is not None:
            refusals = self.flush(raise_on_refusal=False)
            own = [r for r in refusals if r.object_id == id]
            others = [r for r in refusals if r.object_id != id]
            # Keep the refusals of other objects for the next flush.
            self._refusals = others + self._refusals
            if own:
                raise ViewerRefusedError(own)
            self.apply_colour(id, colour)
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
        self._revisions[id] = revision
        self._inventory[id] = InventoryEntry(id, "Group", revision)
        return revision

    def remove(self, id: str) -> None:
        """Remove one object from the scene. The revision count of the id is dropped.

        Parameters
        ----------
        id : str
            The identifier of the object.
        """
        self._ensure_connected()
        self._send_frame(encode.remove_frame(id))
        self._inventory.pop(id, None)
        self._revisions.pop(id, None)
        self._inventory_revisions.pop(id, None)

    def clear(self) -> None:
        """Remove all objects from the scene."""
        self._ensure_connected()
        self._send_frame(encode.clear_frame())
        self._inventory.clear()
        self._revisions.clear()
        self._inventory_revisions.clear()

    def inventory(self) -> list[InventoryEntry]:
        """Return the known scene.

        The list starts with the inventory from the hello frame. ``send`` and
        ``send_group`` add or update an entry, ``remove`` drops one, and ``clear``
        drops all. A new connection replaces the list with the hello list.

        Returns
        -------
        list of InventoryEntry
            The entries.
        """
        return list(self._inventory.values())

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

    # -- selection ------------------------------------------------------------

    def select(
        self,
        object_id: str,
        *,
        cells: Sequence[int] | np.ndarray | None = None,
        vertices: Sequence[int] | np.ndarray | None = None,
        revision: int | None = None,
    ) -> dict[str, int]:
        """Select cells and vertices of one object in the viewer.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        cells : sequence of int, optional
            The cell indices. Each index is from 0 to 2**32-1.
        vertices : sequence of int, optional
            The vertex indices. Each index is from 0 to 2**32-1.
        revision : int, optional
            The revision that the selection is for. The viewer refuses an old revision.

        Returns
        -------
        dict of str to int
            The number of selected cells (key ``cell``) and vertices (key ``vertex``).

        Raises
        ------
        loop_link.errors.InvalidInputError
            The id is not valid, no index list is given, an index list is not valid,
            or the viewer refused the indices.
        loop_link.errors.StaleRevisionError
            The revision is old.
        loop_link.errors.MissingFeatureError
            The viewer does not support the selection.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        if not isinstance(object_id, str) or not object_id:
            raise InvalidInputError("The object_id must be a string that is not empty.")
        if cells is None and vertices is None:
            raise InvalidInputError(
                "Give cells or vertices. Use an empty list to clear the selection."
            )
        cell_array = None if cells is None else _indices(cells, "cells")
        vertex_array = None if vertices is None else _indices(vertices, "vertices")
        tool = self._tool()
        frame = tool.request(
            lambda request: encode.selection_set_frame(
                request, object_id, cell_array, vertex_array, revision
            ),
            until=("selection.result",),
            feature=_SELECTION_FEATURE,
            object_id=object_id,
        )
        counts = {"cell": 0, "vertex": 0}
        for entry in frame.header.get("entries") or []:
            name = _KIND_NAMES.get(entry.get("kind"))
            if name is not None:
                counts[name] = int(entry.get("count", 0))
        return counts

    def selection(self) -> Selection:
        """Read the selection in the viewer.

        Returns
        -------
        Selection
            The selection. It is empty when nothing is selected.

        Raises
        ------
        loop_link.errors.MissingFeatureError
            The viewer does not support the selection.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        tool = self._tool()
        frame = tool.request(
            encode.selection_get_frame,
            until=("selection.result",),
            feature=_SELECTION_FEATURE,
        )
        return parse_selection_result(frame)

    def watch_selection(
        self,
        callback: Callable[[Selection], None],
        *,
        interval: float = 0.5,
        initial: bool = False,
    ) -> SelectionWatcher:
        """Call a function each time the selection in the viewer changes.

        Parameters
        ----------
        callback : callable
            The function to call with the new ``Selection``.
        interval : float
            The time in seconds between two polls. It must be more than zero.
        initial : bool
            If true, call the function once with the current selection when it is
            not empty.

        Returns
        -------
        SelectionWatcher
            The running watcher. Call its ``stop`` method to end it.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The interval is not more than zero.
        loop_link.errors.MissingFeatureError
            The viewer does not support the selection.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        if not interval > 0:
            raise InvalidInputError("The interval must be more than zero.")
        last = self.selection()
        if initial and not last.empty:
            callback(last)
        watcher = SelectionWatcher(weakref.ref(self), callback, interval, last)
        self._watchers.append(watcher)
        watcher.start()
        return watcher

    # -- filters --------------------------------------------------------------

    def filter_range(
        self,
        object_id: str,
        array: str,
        minimum: float,
        maximum: float,
        *,
        keep_nan: bool = False,
        wait: bool = True,
        cancel: threading.Event | None = None,
        timeout: float = 60.0,
        progress: Callable[[float], None] | None = None,
    ) -> FilterResult:
        """Keep the cells or vertices with an array value in a range.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        array : str
            The name of a number array of the object.
        minimum : float
            The low limit of the range.
        maximum : float
            The high limit of the range. It is not less than ``minimum``.
        keep_nan : bool
            If true, keep the items that have no value (NaN).
        wait : bool
            If true, wait for the end of the filter build. If false, send the last edit and
            return an empty result at once. A refusal of that edit is not seen.
        cancel : threading.Event, optional
            When it is set, the viewer cancels the build and the call raises ``CancelledError``.
        timeout : float
            The time in seconds to wait for the next frame of the viewer. It must be more
            than zero.
        progress : callable, optional
            Called with the fraction of each progress frame of the build.

        Returns
        -------
        FilterResult
            The counts of the build. A field is ``None`` when the viewer gives no value.

        Raises
        ------
        loop_link.errors.CancelledError
            The ``cancel`` event was set.
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the filter.
        loop_link.errors.MissingFeatureError
            The viewer does not support the filter commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        array = check_text(array, "array")
        low = check_number(minimum, "minimum")
        high = check_number(maximum, "maximum")
        if low > high:
            raise InvalidInputError("'minimum' must not be more than 'maximum'.")
        keep_nan = _check_bool(keep_nan, "keep_nan")
        _check_wait(wait, timeout)
        tool = self._tool()
        self._filter_add(tool, object_id, array, cancel, timeout)
        return _filter_edits(
            tool,
            object_id,
            (
                {
                    "kind": "set_range",
                    "array": array,
                    "min": format_number(low),
                    "max": format_number(high),
                    "keep_nan": keep_nan,
                },
                {"kind": "set_on", "array": array, "on": True},
            ),
            wait=wait,
            cancel=cancel,
            timeout=timeout,
            progress=progress,
        )

    def filter_values(
        self,
        object_id: str,
        array: str,
        values: Sequence[int] | np.ndarray,
        *,
        exclude: bool = False,
        wait: bool = True,
        cancel: threading.Event | None = None,
        timeout: float = 60.0,
        progress: Callable[[float], None] | None = None,
    ) -> FilterResult:
        """Keep the cells or vertices with an array value in a list of codes.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        array : str
            The name of a category array of the object.
        values : sequence of int
            The category codes. Each code is from 0 to 2**32-1.
        exclude : bool
            If true, remove the items with these codes and keep the others.
        wait : bool
            If true, wait for the end of the filter build. If false, send the last edit and
            return an empty result at once. A refusal of that edit is not seen.
        cancel : threading.Event, optional
            When it is set, the viewer cancels the build and the call raises ``CancelledError``.
        timeout : float
            The time in seconds to wait for the next frame of the viewer. It must be more
            than zero.
        progress : callable, optional
            Called with the fraction of each progress frame of the build.

        Returns
        -------
        FilterResult
            The counts of the build. A field is ``None`` when the viewer gives no value.

        Raises
        ------
        loop_link.errors.CancelledError
            The ``cancel`` event was set.
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the filter.
        loop_link.errors.MissingFeatureError
            The viewer does not support the filter commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        array = check_text(array, "array")
        codes = check_codes(values, "values")
        exclude = _check_bool(exclude, "exclude")
        _check_wait(wait, timeout)
        tool = self._tool()
        self._filter_add(tool, object_id, array, cancel, timeout)
        return _filter_edits(
            tool,
            object_id,
            (
                {"kind": "set_values", "array": array, "values": codes},
                {"kind": "set_exclude", "array": array, "exclude": exclude},
                {"kind": "set_on", "array": array, "on": True},
            ),
            wait=wait,
            cancel=cancel,
            timeout=timeout,
            progress=progress,
        )

    def filter_enable(
        self,
        object_id: str,
        array: str,
        on: bool = True,
        *,
        wait: bool = True,
        cancel: threading.Event | None = None,
        timeout: float = 60.0,
        progress: Callable[[float], None] | None = None,
    ) -> FilterResult:
        """Turn the filter of one array on or off.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        array : str
            The name of the array of the filter.
        on : bool
            True to turn the filter on, false to turn it off.
        wait : bool
            If true, wait for the end of the filter build. If false, send the edit and
            return an empty result at once.
        cancel : threading.Event, optional
            When it is set, the viewer cancels the build and the call raises ``CancelledError``.
        timeout : float
            The time in seconds to wait for the next frame of the viewer. It must be more
            than zero.
        progress : callable, optional
            Called with the fraction of each progress frame of the build.

        Returns
        -------
        FilterResult
            The counts of the build. A field is ``None`` when the viewer gives no value.

        Raises
        ------
        loop_link.errors.CancelledError
            The ``cancel`` event was set.
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the filter commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        array = check_text(array, "array")
        on = _check_bool(on, "on")
        _check_wait(wait, timeout)
        return _filter_edits(
            self._tool(),
            object_id,
            ({"kind": "set_on", "array": array, "on": on},),
            wait=wait,
            cancel=cancel,
            timeout=timeout,
            progress=progress,
        )

    def filter_remove(
        self,
        object_id: str,
        array: str,
        *,
        wait: bool = True,
        cancel: threading.Event | None = None,
        timeout: float = 60.0,
        progress: Callable[[float], None] | None = None,
    ) -> FilterResult:
        """Remove the filter of one array.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        array : str
            The name of the array of the filter.
        wait : bool
            If true, wait for the end of the filter build. If false, send the edit and
            return an empty result at once.
        cancel : threading.Event, optional
            When it is set, the viewer cancels the build and the call raises ``CancelledError``.
        timeout : float
            The time in seconds to wait for the next frame of the viewer. It must be more
            than zero.
        progress : callable, optional
            Called with the fraction of each progress frame of the build.

        Returns
        -------
        FilterResult
            The counts of the build. A field is ``None`` when the viewer gives no value.

        Raises
        ------
        loop_link.errors.CancelledError
            The ``cancel`` event was set.
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the filter commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        array = check_text(array, "array")
        _check_wait(wait, timeout)
        return _filter_edits(
            self._tool(),
            object_id,
            ({"kind": "remove", "array": array},),
            wait=wait,
            cancel=cancel,
            timeout=timeout,
            progress=progress,
        )

    def cancel_filter(self, object_id: str) -> None:
        """Cancel the filter build that runs for one object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The id is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the filter commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        _filter_edit(self._tool(), object_id, {"kind": "cancel"})

    @staticmethod
    def _filter_add(
        tool: ToolConnection,
        object_id: str,
        array: str,
        cancel: threading.Event | None,
        timeout: float,
    ) -> None:
        """Add a filter. A refusal is not an error: the filter can exist already."""
        try:
            _filter_edit(tool, object_id, {"kind": "add", "array": array}, cancel, timeout)
        except InvalidInputError:
            pass

    # -- colour ---------------------------------------------------------------

    def colour_map_names(self) -> tuple[str, ...]:
        """Return the names of the colour maps of the viewer.

        The names are kept until the viewer is closed.

        Returns
        -------
        tuple of str
            The names.

        Raises
        ------
        loop_link.errors.MissingFeatureError
            The viewer does not give the names.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        names = self._read_colour_map_names(self._tool())
        if names is None:
            raise MissingFeatureError("colour map names (geoviewer protocol 1.18)")
        return names

    def set_colour_property(self, object_id: str, array: str | None) -> None:
        """Set the array that colours an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        array : str or None
            The name of the array. ``None`` shows no array.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        body: dict[str, Any] = {"id": object_id}
        if array is not None:
            body["array"] = check_text(array, "array")
        run_command(self._tool(), "colour.set_active", body, object_id=object_id)

    def set_colour_map(self, object_id: str, name: str) -> None:
        """Set the colour map of an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        name : str
            The name of the colour map.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, the name is not known, or the viewer refused
            the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        name = check_text(name, "name")
        tool = self._tool()
        self._check_colour_map(tool, name)
        run_command(tool, "colour.set_map", {"id": object_id, "name": name}, object_id=object_id)

    def set_colour_range(self, object_id: str, minimum: float, maximum: float) -> None:
        """Set the colour range of an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        minimum : float
            The low limit of the range.
        maximum : float
            The high limit of the range. It is more than ``minimum``.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        low = check_number(minimum, "minimum")
        high = check_number(maximum, "maximum")
        if not low < high:
            raise InvalidInputError("'minimum' must be less than 'maximum'.")
        body = {"id": object_id, "min": format_number(low), "max": format_number(high)}
        run_command(self._tool(), "colour.set_range", body, object_id=object_id)

    def clear_colour_range(self, object_id: str) -> None:
        """Remove the colour range of an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The id is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        run_command(self._tool(), "colour.clear_range", {"id": object_id}, object_id=object_id)

    def set_colour_scale(self, object_id: str, scale: str) -> None:
        """Set the colour scale of an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        scale : str
            The scale, ``"linear"`` or ``"log"``.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour scale.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        if scale not in _COLOUR_SCALES:
            raise InvalidInputError(f"'scale' must be 'linear' or 'log'. Got {scale!r}.")
        body = {"object": object_id, "scale": scale}
        run_command(self._tool(), "colour.set_scale", body, object_id=object_id)

    def set_object_colour(self, object_id: str, colour: Sequence[int] | None) -> None:
        """Set the plain colour of an object.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        colour : sequence of int or None
            The red, green and blue values from 0 to 255. ``None`` removes the plain colour.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused the command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        value = None if colour is None else list(check_rgb(colour, "colour"))
        body = {"id": object_id, "colour": value}
        run_command(self._tool(), "object.set_colour", body, object_id=object_id)

    def apply_colour(self, object_id: str, settings: ColourSettings) -> None:
        """Apply the colour settings of an object.

        The commands go in this order: array, colour map, range, scale, plain colour.
        A field that is ``None`` sends no command.

        Parameters
        ----------
        object_id : str
            The identifier of the object.
        settings : loop_link.commands.ColourSettings
            The colour settings.

        Raises
        ------
        loop_link.errors.InvalidInputError
            An argument is not valid, or the viewer refused a command.
        loop_link.errors.MissingFeatureError
            The viewer does not support the colour commands.
        loop_link.errors.ViewerClosedError
            The connection is closed.
        """
        object_id = check_text(object_id, "object_id")
        if not isinstance(settings, ColourSettings):
            raise InvalidInputError("The settings must be a ColourSettings object.")
        if settings.array is not None:
            self.set_colour_property(object_id, settings.array)
        if settings.colour_map is not None:
            self.set_colour_map(object_id, settings.colour_map)
        if settings.range is not None:
            self.set_colour_range(object_id, *settings.range)
        if settings.scale is not None:
            self.set_colour_scale(object_id, settings.scale)
        if settings.plain is not None:
            self.set_object_colour(object_id, settings.plain)

    def _read_colour_map_names(self, tool: ToolConnection) -> tuple[str, ...] | None:
        """Return the colour map names. Ask the viewer at the first call only."""
        if self._colour_map_names is None:
            self._colour_map_names = list_colour_map_names(tool)
        return self._colour_map_names

    def _check_colour_map(self, tool: ToolConnection, name: str) -> None:
        """Raise when the viewer has names and ``name`` is not one of them."""
        names = self._read_colour_map_names(tool)
        if names is not None and name not in names:
            raise InvalidInputError(
                f"The viewer has no colour map '{name}'. The names are: {', '.join(names)}."
            )

    # -- fetching -------------------------------------------------------------

    def fetch(
        self,
        id: str,
        *,
        progress: Callable[[float], None] | None = None,
        cancel: threading.Event | None = None,
    ) -> FetchResult:
        """Read one object from the viewer.

        This method uses the tool connection. It does not open the data connection.

        Parameters
        ----------
        id : str
            The identifier of the object.
        progress : callable, optional
            The function to call with the fraction done, from 0 to 1.
        cancel : threading.Event, optional
            When this event is set, the fetch stops.

        Returns
        -------
        FetchResult
            The object, with its id, type, revision and name.

        Raises
        ------
        loop_link.errors.InvalidInputError
            The id is not valid, or the viewer does not know the id.
        loop_link.errors.CancelledError
            The ``cancel`` event was set.
        loop_link.errors.MissingFeatureError
            The viewer does not support the fetch.
        loop_link.errors.ViewerClosedError
            The connection is closed, or the viewer sent an object that is not valid.
        loop_link.errors.UnsupportedObjectError
            The viewer sent an object type that this package cannot read.
        loop_link.errors.NoViewerError
            No viewer was found.
        loop_link.errors.TokenError
            The viewer refused the token.
        loop_link.errors.ObjectTooLargeError
            The object is larger than the limit.
        """
        if not isinstance(id, str) or not id:
            raise InvalidInputError("The id must be a string that is not empty.")
        tool = self._tool()
        return fetch_object(tool, id, progress=progress, cancel=cancel)

    # -- closing --------------------------------------------------------------

    def release_data(self) -> None:
        """Close the connection so another client can own the scene.

        A later ``send`` opens a new connection. The revision numbers stay. The
        tool connection stays open. A later ``send`` also opens a new connection
        after the connection fails.
        """
        self._stop_thread()
        self._drop_client()

    def close(self) -> None:
        """Close the connection and stop all selection watchers. The revision numbers stay."""
        watchers, self._watchers = self._watchers, []
        for watcher in watchers:
            watcher.stop()
        with self._tool_lock:
            self._closed = True
            tool, self._tool_connection = self._tool_connection, None
            self._colour_map_names = None
        if tool is not None:
            tool.close()
        self._stop_thread()
        self._drop_client()
        self._dead = True
        with self._address_lock:
            self._address = None

    # -- helpers --------------------------------------------------------------

    def _tool(self) -> ToolConnection:
        """Return the tool connection. Make it at the first call."""
        if self._closed:
            raise ViewerClosedError
        with self._tool_lock:
            if self._closed:
                raise ViewerClosedError
            if self._tool_connection is None:
                self._tool_connection = ToolConnection(
                    timeout=self._timeout, resolve=self._resolve_address
                )
            return self._tool_connection

    def _resolve_address(self) -> discovery.ViewerAddress:
        """Return the address of the viewer.

        While the data connection or the tool connection is open, the known address
        is used again, so that both connections go to the same viewer. Otherwise the
        viewer is found again.
        """
        with self._address_lock:
            tool = self._tool_connection
            in_use = (self._client is not None and not self._dead) or (
                tool is not None and tool.connected
            )
            if self._address is None or not in_use:
                self._address = discovery.find_viewer(
                    port=self._port,
                    token=self._token,
                    host=self._host,
                    settings=self._settings_dir,
                )
            return self._address

    def _forget_address(self) -> None:
        """Forget the known address and close the tool connection that uses it."""
        with self._tool_lock:
            tool = self._tool_connection
            self._colour_map_names = None
        if tool is not None:
            tool.close()
        with self._address_lock:
            self._address = None

    def _ensure_connected(self) -> None:
        if self._closed:
            raise ViewerClosedError
        if self._client is None:
            self.connect()
        elif self._dead or self._client.closed:
            # The last connection failed (a client closes itself after a broken send): open a new one. If the viewer is gone, give
            # the error of a closed viewer, as the failed send did.
            try:
                self.connect()
            except NoViewerError as error:
                raise ViewerClosedError(
                    "The viewer closed the connection and no viewer was found again. "
                    "Start geoviewer and turn the link on, then send again."
                ) from error

    def _drop_client(self) -> None:
        """Close the socket, if there is one."""
        with self._lock:
            client, self._client = self._client, None
            self._shared["client"] = None
        if client is not None:
            _close_client(client)

    def _next_revision(self, id: str) -> int:
        return max(self._revisions.get(id, 0), self._inventory_revisions.get(id, 0)) + 1

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
            try:
                header = protocol.decode_frame(data).header
            except ProtocolError as error:
                self._dead = True
                _close_client(client)
                raise ViewerClosedError(
                    "The viewer sent a frame that is not valid. The connection is closed."
                ) from error
            if header.get("type") == "error" and header.get("object") is not None:
                self._refusals.append(
                    Refusal(
                        str(header["object"]),
                        str(header.get("code", "")),
                        str(header.get("message", "")),
                    )
                )
            wait = 0.001

    def _heartbeat_tick(self) -> bool:
        """Do one keepalive step. Return False when the thread must stop."""
        if not self._lock.acquire(blocking=False):
            return True
        try:
            client = self._client
            if client is None or self._dead:
                return False
            try:
                self._drain_locked(client, 0.001)
            except _LINK_ERRORS + (ViewerClosedError,):
                self._dead = True
                return False
            if time.monotonic() - self._last_send < self._heartbeat_interval:
                return True
            self._seq += 1
            try:
                client.send_binary(encode.heartbeat_frame(self._seq))
                self._last_send = time.monotonic()
            except _LINK_ERRORS:
                self._dead = True
                return False
            return True
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


def _indices(values: Any, name: str) -> np.ndarray:  # noqa: ANN401
    """Return a 1-D ``uint32`` array from index values. Raise InvalidInputError if bad."""
    try:
        array = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise InvalidInputError(f"The {name} must be a list of integers.") from error
    if array.size == 0 and array.ndim == 1:
        return np.empty(0, dtype=np.uint32)
    if array.dtype.kind not in "iu":
        raise InvalidInputError(f"The {name} must be integers. Do not use float or bool values.")
    if array.ndim != 1:
        raise InvalidInputError(f"The {name} must be a 1-D list of integers.")
    if array.min() < 0 or array.max() > _INDEX_MAX:
        raise InvalidInputError(f"The {name} must be from 0 to {_INDEX_MAX}.")
    return array.astype(np.uint32)


def _check_bool(value: object, name: str) -> bool:
    """Return a value that must be a ``bool``. Raise InvalidInputError if it is not."""
    if not isinstance(value, bool):
        raise InvalidInputError(f"'{name}' must be True or False. Got {value!r}.")
    return value


def _check_wait(wait: object, timeout: object) -> None:
    """Check the ``wait`` and ``timeout`` arguments of a filter call."""
    _check_bool(wait, "wait")
    if check_number(timeout, "timeout") <= 0:
        raise InvalidInputError("'timeout' must be more than zero.")


def _filter_edit(
    tool: ToolConnection,
    object_id: str,
    edit: dict[str, Any],
    cancel: threading.Event | None = None,
    timeout: float | None = None,
    progress: Callable[[float], None] | None = None,
) -> Mapping[str, Any]:
    """Send one ``filter.edit`` command and wait for its result."""
    return run_command(
        tool,
        "filter.edit",
        {"object": object_id, "edit": edit},
        object_id=object_id,
        cancel=cancel,
        timeout=timeout,
        progress=progress,
    )


def _filter_edits(
    tool: ToolConnection,
    object_id: str,
    edits: Sequence[dict[str, Any]],
    *,
    wait: bool,
    cancel: threading.Event | None,
    timeout: float,
    progress: Callable[[float], None] | None,
) -> FilterResult:
    """Send the edits of a filter call and return the result of the last one.

    With ``wait`` false, the last edit goes with ``post`` and the result is empty.
    """
    result: Mapping[str, Any] = {}
    last = len(edits) - 1
    for index, edit in enumerate(edits):
        if index == last and not wait:
            body = {"object": object_id, "edit": edit}
            tool.post(lambda number: encode.command_frame(number, "filter.edit", body))
            return FilterResult()
        result = _filter_edit(tool, object_id, edit, cancel, timeout, progress)
    return FilterResult.from_result(result)


def _close_client(client: WebSocketClient) -> None:
    """Close a client. An error on close does not matter."""
    try:
        client.close()
    except Exception:
        pass


def _raise_hello_error(header: dict[str, Any]) -> None:
    """Raise the right error for an ``error`` frame during the hello."""
    code = header.get("code")
    message = str(header.get("message", ""))
    if code == "busy":
        raise BusyError(message or "The viewer is busy. Wait a moment and try again.")
    if code == "unsupported_version":
        supported = header.get("supported")
        feature = "protocol version"
        if supported is not None:
            feature = f"protocol 1.x (the viewer supports: {supported})"
        raise MissingFeatureError(feature)
    raise ViewerRefusedError([Refusal(str(header.get("object", "")), str(code), message)])


def _parse_hello(hello: dict[str, Any]) -> tuple[tuple[float, float, float] | None, str]:
    """Return the origin and the CRS of a ``hello`` frame. Raise ProtocolError if it is bad."""
    try:
        origin = hello.get("origin")
        crs = hello.get("crs")
        return (
            _triple(origin) if origin is not None else None,
            "" if crs is None else str(crs),
        )
    except (TypeError, ValueError, AttributeError) as error:
        raise ProtocolError("The viewer sent a hello frame that is not valid.") from error


def _parse_inventory(header: dict[str, Any]) -> list[InventoryEntry]:
    """Return the entries of a ``scene.inventory`` frame. Raise ProtocolError if it is bad."""
    try:
        return [
            InventoryEntry(
                id=str(item["id"]),
                object_type=str(item.get("object_type", "")),
                revision=int(item.get("revision", 1)),
            )
            for item in header.get("entries", [])
        ]
    except (TypeError, ValueError, KeyError, AttributeError) as error:
        raise ProtocolError("The viewer sent an inventory frame that is not valid.") from error


def _read_early(client: WebSocketClient, wait: float = 0.05) -> None:
    """Read frames that wait already. Raise the hello error if one is an ``error`` frame."""
    for _ in range(4):
        data = client.receive(wait)
        if data is None:
            return
        header = protocol.decode_frame(data).header
        if header.get("type") == "error":
            _raise_hello_error(header)


def _heartbeat_loop(
    ref: weakref.ReferenceType[Viewer], stop: threading.Event, interval: float
) -> None:
    """Run the keepalive. The thread holds only a weak reference to the viewer."""
    while not stop.wait(interval):
        viewer = ref()
        if viewer is None or not viewer._heartbeat_tick():
            return
        del viewer


def _release_shared(shared: dict[str, Any]) -> None:
    """Stop the thread and close the socket of a viewer that is gone."""
    shared["stop"].set()
    client = shared.get("client")
    if client is not None:
        _close_client(client)
