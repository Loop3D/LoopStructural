"""The tool connection: one request and its reply frames on the ``/tool`` path of the viewer."""

from __future__ import annotations

import threading
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, NoReturn

from loop_link import discovery, encode, protocol
from loop_link._websocket import WebSocketClient, WebSocketClosed, WebSocketError
from loop_link.errors import (
    BusyError,
    CancelledError,
    InvalidInputError,
    MissingFeatureError,
    Refusal,
    StaleRevisionError,
    ViewerClosedError,
    ViewerRefusedError,
)
from loop_link.protocol import Frame, ProtocolError

_LINK_ERRORS = (OSError, WebSocketClosed, WebSocketError)
_POLL = 0.1


def raise_for_request_error(error: Mapping[str, Any], *, object_id: str | None) -> NoReturn:
    """Raise the loop_link error that matches a ``request.error`` body.

    Parameters
    ----------
    error : mapping
        The ``error`` body of the frame, with the keys ``code``, ``message`` and ``field``.
    object_id : str, optional
        The identifier of the object that the request is about.

    Raises
    ------
    loop_link.errors.StaleRevisionError
        The code is ``invalid_input`` and the field is ``revision``.
    loop_link.errors.InvalidInputError
        The code is ``invalid_input`` and the field is another field.
    loop_link.errors.MissingFeatureError
        The code is ``needs_window``.
    loop_link.errors.CancelledError
        The code is ``cancelled``.
    loop_link.errors.BusyError
        The code is ``busy``.
    loop_link.errors.ViewerRefusedError
        The code is another code.
    """
    code = str(error.get("code", ""))
    message = str(error.get("message", ""))
    if code == "invalid_input":
        if error.get("field") == "revision":
            raise StaleRevisionError(object_id or "")
        raise InvalidInputError(f"The viewer refused the request: {message}")
    if code == "needs_window":
        raise MissingFeatureError("a viewer window (this viewer has no window)")
    if code == "cancelled":
        raise CancelledError
    if code == "busy":
        raise BusyError
    raise ViewerRefusedError([Refusal(object_id or "", code, message)])


class ToolConnection:
    """A connection to the ``/tool`` path of the viewer.

    The connection sends one request at a time. It opens when the first request needs it.

    Parameters
    ----------
    port : int, optional
        The port of the viewer. The default is to search for it.
    token : str, optional
        The access token. The default is to read it from the environment or the settings.
    host : str
        The host of the viewer.
    settings_dir : pathlib.Path, optional
        The settings folder of the viewer.
    timeout : float
        The time in seconds to wait for the connection and for each reply frame.
    """

    def __init__(
        self,
        *,
        port: int | None = None,
        token: str | None = None,
        host: str = "127.0.0.1",
        settings_dir: Path | None = None,
        timeout: float = 5.0,
    ) -> None:
        self._port = port
        self._token = token
        self._host = host
        self._settings_dir = settings_dir
        self._timeout = timeout
        self._lock = threading.Lock()
        self._client: WebSocketClient | None = None
        self._number = 0

    @property
    def connected(self) -> bool:
        """True while the connection is open."""
        return self._client is not None

    def connect(self) -> ToolConnection:
        """Connect to the ``/tool`` path of the viewer. No ``hello`` is sent.

        Raises
        ------
        loop_link.errors.NoViewerError
            No viewer was found.
        loop_link.errors.TokenError
            The viewer refused the token.

        Returns
        -------
        ToolConnection
            The connection itself.
        """
        with self._lock:
            self._connect_locked()
        return self

    def close(self) -> None:
        """Close the connection. It is safe to call it more than one time."""
        with self._lock:
            self._drop_locked()

    def request(
        self,
        build: Callable[[int], bytes],
        *,
        until: Sequence[str],
        feature: str,
        object_id: str | None = None,
        on_frame: Callable[[Frame], None] | None = None,
        progress: Callable[[float], None] | None = None,
        cancel: threading.Event | None = None,
        timeout: float | None = None,
    ) -> Frame:
        """Send one request and read frames until the reply ends.

        Parameters
        ----------
        build : callable
            Gets the request number and returns the frame to send.
        until : sequence of str
            The frame types that end the request. The matching frame is returned.
        feature : str
            The name of the feature, for the error when the viewer does not know the message.
        object_id : str, optional
            The identifier of the object that the request is about, for error messages.
        on_frame : callable, optional
            Called with each frame of the request, except ``command.progress`` and
            ``request.error`` frames.
        progress : callable, optional
            Called with the fraction of each ``command.progress`` frame.
        cancel : threading.Event, optional
            When it is set, one ``command.cancel`` frame is sent.
        timeout : float, optional
            The time in seconds to wait for the next frame. The default is the
            timeout of the connection.

        Returns
        -------
        loop_link.protocol.Frame
            The frame whose type is in ``until``.

        Raises
        ------
        loop_link.errors.NoViewerError
            No viewer was found.
        loop_link.errors.TokenError
            The viewer refused the token.
        loop_link.errors.ViewerClosedError
            The connection closed, the viewer did not answer in time, or it sent a bad frame.
        loop_link.errors.MissingFeatureError
            The viewer does not know the message or the protocol version.
        loop_link.errors.CancelledError
            The request was cancelled.
        loop_link.errors.ViewerRefusedError
            The viewer refused the request.
        """
        wait = self._timeout if timeout is None else timeout
        with self._lock:
            client = self._connect_locked()
            self._number += 1
            number = self._number
            try:
                client.send_binary(build(number))
                return self._read_locked(
                    client, number, until, feature, object_id, on_frame, progress, cancel, wait
                )
            except _LINK_ERRORS as error:
                self._drop_locked()
                raise ViewerClosedError from error

    # -- helpers --------------------------------------------------------------

    def _connect_locked(self) -> WebSocketClient:
        """Return the open client, or open one. The caller holds the lock."""
        if self._client is not None:
            return self._client
        address = discovery.find_viewer(
            port=self._port,
            token=self._token,
            host=self._host,
            settings=self._settings_dir,
        )
        self._client = WebSocketClient.connect(
            address.host, address.port, token=address.token, timeout=self._timeout, path="/tool"
        )
        return self._client

    def _drop_locked(self) -> None:
        """Close the socket, if there is one. The caller holds the lock."""
        client, self._client = self._client, None
        if client is not None:
            try:
                client.close()
            except Exception:  # noqa: BLE001, S110
                pass

    def _read_locked(
        self,
        client: WebSocketClient,
        number: int,
        until: Sequence[str],
        feature: str,
        object_id: str | None,
        on_frame: Callable[[Frame], None] | None,
        progress: Callable[[float], None] | None,
        cancel: threading.Event | None,
        wait: float,
    ) -> Frame:
        """Read the frames of one request. The caller holds the lock."""
        last_frame = time.monotonic()
        cancel_sent = False
        while True:
            if cancel is not None and cancel.is_set() and not cancel_sent:
                client.send_binary(encode.cancel_frame(number))
                cancel_sent = True
            data = client.receive(_POLL)
            if data is None:
                if time.monotonic() - last_frame >= wait:
                    self._drop_locked()
                    raise ViewerClosedError(
                        "The viewer did not answer in time. Check the viewer and retry."
                    )
                continue
            try:
                frame = protocol.decode_frame(data)
            except ProtocolError as error:
                self._drop_locked()
                raise ViewerClosedError(
                    "The viewer sent a frame that is not valid. The connection is closed."
                ) from error
            header = frame.header
            kind = header.get("type")
            if kind == "error":
                self._raise_error_frame(header, feature)
            if header.get("request") != number:
                continue
            last_frame = time.monotonic()
            if kind == "request.error":
                raise_for_request_error(header.get("error") or {}, object_id=object_id)
            if kind == "command.progress":
                fraction = header.get("fraction")
                if progress is not None and isinstance(fraction, (int, float)):
                    progress(float(fraction))
                continue
            if on_frame is not None:
                on_frame(frame)
            if kind in until:
                return frame

    @staticmethod
    def _raise_error_frame(header: Mapping[str, Any], feature: str) -> NoReturn:
        """Raise the error for an ``error`` frame, which has no request number."""
        code = str(header.get("code", ""))
        if code == "unknown_message":
            raise MissingFeatureError(feature)
        if code == "unsupported_version":
            raise MissingFeatureError("protocol " + str(header.get("supported")))
        raise ViewerRefusedError([Refusal("", code, str(header.get("message", "")))])
