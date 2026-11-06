"""A small WebSocket client for the link to the viewer.

It uses the standard library only, because some host applications do not ship a
WebSocket package for every platform. It does what the link needs: one
connection with an ``Authorization`` header, binary messages in both
directions, and the control frames (ping, pong, close). It does not do TLS.
"""

from __future__ import annotations

import base64
import hashlib
import os
import socket
import struct
import time

from loop_link.errors import LoopLinkError, NoViewerError, TokenError
from loop_link.protocol import MAX_FRAME_BYTES, MAX_HEADER_BYTES

_GUID = b"258EAFA5-E914-47DA-95CA-C5AB0DC85B11"

OP_CONTINUATION = 0x0
OP_TEXT = 0x1
OP_BINARY = 0x2
OP_CLOSE = 0x8
OP_PING = 0x9
OP_PONG = 0xA

# The largest incoming message: one protocol frame, which is a header and a body.
MAX_INCOMING_BYTES: int = MAX_FRAME_BYTES + MAX_HEADER_BYTES
_HANDSHAKE_MAX_BYTES = 16 * 1024
# Time for the rest of a frame after its first bytes arrived
_FRAME_TIMEOUT = 10.0
# Time to send one frame. A receive can leave a short timeout on the socket,
# so each send sets this first.
_SEND_TIMEOUT = 60.0
# The smallest socket timeout: a timeout of zero makes the socket non-blocking.
_MIN_WAIT = 1e-3


class WebSocketError(LoopLinkError):
    """The connection failed, or the server broke the protocol."""


class WebSocketClosed(WebSocketError):
    """The server closed the connection."""


def _mask(data: bytes, key: bytes) -> bytes:
    """Apply the WebSocket mask to data.

    Parameters
    ----------
    data : bytes
        The bytes to mask or unmask.
    key : bytes
        The 4-byte mask key.

    Returns
    -------
    bytes
        The data after the XOR with the repeated key.
    """
    if not data:
        return data
    repeated = (key * (len(data) // 4 + 1))[: len(data)]
    value = int.from_bytes(data, "little") ^ int.from_bytes(repeated, "little")
    return value.to_bytes(len(data), "little")


class WebSocketClient:
    """One client connection. Not thread safe: use it from one thread.

    Parameters
    ----------
    sock : socket.socket
        A connected socket. Use :meth:`connect` to make a client.

    Attributes
    ----------
    closed : bool
        True after the connection is closed.
    """

    def __init__(self, sock: socket.socket) -> None:
        self._sock = sock
        self._buffer = b""
        self._partial = False
        self._fragmenting = False
        self._fragment_opcode = OP_BINARY
        self._fragments = bytearray()
        self.closed = False

    @classmethod
    def connect(
        cls,
        host: str,
        port: int,
        *,
        token: str = "",
        timeout: float = 5.0,
        path: str = "/",
    ) -> WebSocketClient:
        """Open the connection and do the handshake.

        Parameters
        ----------
        host : str
            The host name or address of the viewer.
        port : int
            The TCP port of the viewer.
        token : str, optional
            The access token. It is sent as ``Authorization: Bearer <token>``.
        timeout : float, optional
            The time in seconds to wait for the TCP connection.
        path : str, optional
            The path of the WebSocket request.

        Returns
        -------
        WebSocketClient
            The connected client.

        Raises
        ------
        NoViewerError
            If the TCP connection is refused or times out.
        TokenError
            If the viewer answers with HTTP 401.
        WebSocketError
            If the handshake fails for another reason.
        """
        try:
            sock = socket.create_connection((host, port), timeout=timeout)
        except OSError as error:
            raise NoViewerError(f"Cannot connect to a viewer at {host}:{port}: {error}") from error
        try:
            sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            client = cls(sock)
            client._handshake(host, port, token, path)
        except BaseException:
            sock.close()
            raise
        return client

    def _handshake(self, host: str, port: int, token: str, path: str) -> None:
        """Send the upgrade request and check the answer.

        Parameters
        ----------
        host : str
            The host name for the ``Host`` header.
        port : int
            The port for the ``Host`` header.
        token : str
            The access token. An empty string sends no ``Authorization`` header.
        path : str
            The path of the request.
        """
        key = base64.b64encode(os.urandom(16))
        lines = [
            f"GET {path} HTTP/1.1",
            f"Host: {host}:{port}",
            "Upgrade: websocket",
            "Connection: Upgrade",
            f"Sec-WebSocket-Key: {key.decode()}",
            "Sec-WebSocket-Version: 13",
        ]
        if token:
            lines.append(f"Authorization: Bearer {token}")
        try:
            self._sock.sendall(("\r\n".join(lines) + "\r\n\r\n").encode("ascii"))
            data = b""
            while b"\r\n\r\n" not in data:
                chunk = self._sock.recv(4096)
                if not chunk:
                    raise WebSocketError("The viewer closed the connection during the handshake.")
                data += chunk
                if len(data) > _HANDSHAKE_MAX_BYTES:
                    raise WebSocketError("The handshake answer is too long.")
        except OSError as error:
            raise WebSocketError(f"The handshake failed: {error}") from error
        head, _, rest = data.partition(b"\r\n\r\n")
        self._buffer = rest
        status_line, *header_lines = head.decode("latin-1").split("\r\n")
        parts = status_line.split(" ", 2)
        status = parts[1] if len(parts) > 1 else ""
        if status == "401":
            raise TokenError("The viewer refused the token. Check the token of the link.")
        if status != "101":
            raise WebSocketError(f"The viewer refused the connection: {status_line}")
        headers = {}
        for line in header_lines:
            name, _, value = line.partition(":")
            headers[name.strip().lower()] = value.strip()
        expected = base64.b64encode(hashlib.sha1(key + _GUID).digest()).decode()
        if headers.get("sec-websocket-accept") != expected:
            raise WebSocketError("The server did not answer as a WebSocket server.")

    # -- sending ------------------------------------------------------------

    def _send_frame(self, opcode: int, payload: bytes) -> None:
        """Send one masked frame.

        Parameters
        ----------
        opcode : int
            The frame opcode.
        payload : bytes
            The frame body, before masking.
        """
        length = len(payload)
        head = bytearray([0x80 | opcode])
        if length < 126:
            head.append(0x80 | length)
        elif length < 1 << 16:
            head.append(0x80 | 126)
            head += struct.pack(">H", length)
        else:
            head.append(0x80 | 127)
            head += struct.pack(">Q", length)
        key = os.urandom(4)
        try:
            self._sock.settimeout(_SEND_TIMEOUT)
            self._sock.sendall(bytes(head) + key + _mask(payload, key))
        except OSError as error:
            # This includes a timeout: part of the frame can be sent, so the stream is broken.
            self._abort()
            raise WebSocketError(f"Cannot send to the viewer: {error}") from error

    def send_binary(self, data: bytes) -> None:
        """Send one binary message.

        Parameters
        ----------
        data : bytes
            The message body.

        Raises
        ------
        WebSocketClosed
            If the connection is already closed.
        WebSocketError
            If the send fails. The connection is then closed.
        """
        self._check_open()
        self._send_frame(OP_BINARY, data)

    # -- receiving ----------------------------------------------------------

    def _read_exactly(self, count: int) -> bytes:
        """Read an exact number of bytes from the buffer and the socket.

        Parameters
        ----------
        count : int
            The number of bytes to read.

        Returns
        -------
        bytes
            The bytes that were read.
        """
        while len(self._buffer) < count:
            chunk = self._sock.recv(65536)
            if not chunk:
                self.closed = True
                raise WebSocketClosed("The viewer closed the connection.")
            self._buffer += chunk
        data, self._buffer = self._buffer[:count], self._buffer[count:]
        return data

    def _read_frame(self, wait: float | None) -> tuple[bool, int, bytes]:
        """Read one frame.

        Parameters
        ----------
        wait : float or None
            The time to wait for the first bytes of the frame.

        Returns
        -------
        tuple of (bool, int, bytes)
            The FIN flag, the opcode and the unmasked body.
        """
        # The wait applies to the start of a frame only. A frame that started
        # is read to the end, so a timeout cannot leave half a frame behind.
        self._sock.settimeout(wait)
        first, second = self._read_exactly(2)
        self._partial = True
        self._sock.settimeout(_FRAME_TIMEOUT)
        fin, opcode = bool(first & 0x80), first & 0x0F
        length = second & 0x7F
        if length == 126:
            (length,) = struct.unpack(">H", self._read_exactly(2))
        elif length == 127:
            (length,) = struct.unpack(">Q", self._read_exactly(8))
        if first & 0x70:
            self._abort()
            raise WebSocketError("The viewer sent a frame with a reserved bit set.")
        if opcode >= OP_CLOSE and length > 125:
            self._abort()
            raise WebSocketError("The viewer sent a control frame that is too long.")
        if length > MAX_INCOMING_BYTES:
            self._abort()
            raise WebSocketError("The viewer sent a message that is too large.")
        key = self._read_exactly(4) if second & 0x80 else None
        payload = self._read_exactly(length)
        return fin, opcode, _mask(payload, key) if key else payload

    def receive(self, timeout: float) -> bytes | None:
        """Return the next binary message.

        Ping frames get a pong. A text message is skipped. A binary message in
        parts is joined.

        Parameters
        ----------
        timeout : float
            The time in seconds to wait for a message to start.

        Returns
        -------
        bytes or None
            The message body, or None if nothing arrives in ``timeout`` seconds.
            The socket stays usable after a timeout.

        Raises
        ------
        WebSocketClosed
            If the viewer closes the connection or sends a close frame.
        WebSocketError
            If a message is too large, the stream breaks, or a read fails.
        """
        self._check_open()
        deadline = time.monotonic() + timeout
        try:
            while True:
                self._partial = False
                wait = max(deadline - time.monotonic(), _MIN_WAIT)
                fin, opcode, payload = self._read_frame(wait)
                if opcode == OP_PING:
                    self._send_frame(OP_PONG, payload)
                elif opcode == OP_PONG:
                    continue
                elif opcode == OP_CLOSE:
                    self.closed = True
                    try:
                        self._send_frame(OP_CLOSE, payload[:2])
                    except WebSocketError:
                        pass
                    self._close_socket()
                    raise WebSocketClosed("The viewer closed the connection.")
                else:
                    if opcode == OP_CONTINUATION:
                        if not self._fragmenting:
                            self._abort()
                            raise WebSocketError("The viewer sent an unexpected continuation.")
                    else:
                        if self._fragmenting:
                            self._abort()
                            raise WebSocketError("The viewer sent a new message too early.")
                        self._fragment_opcode = opcode
                        self._fragments = bytearray()
                    self._fragments += payload
                    if len(self._fragments) > MAX_INCOMING_BYTES:
                        self._abort()
                        raise WebSocketError("The viewer sent a message that is too large.")
                    self._fragmenting = not fin
                    if fin:
                        message = bytes(self._fragments)
                        self._fragments = bytearray()
                        if self._fragment_opcode == OP_BINARY:
                            return message
        except (TimeoutError, BlockingIOError) as error:
            if not self._partial:
                # The state of a message in parts stays: a later call goes on with it.
                return None
            # The bytes of the frame that were read are gone, so the next read
            # would start in the middle of a frame.
            self._abort()
            raise WebSocketError("The viewer stopped in the middle of a message.") from error
        except OSError as error:
            self._abort()
            raise WebSocketError(f"Cannot read from the viewer: {error}") from error

    def _check_open(self) -> None:
        """Raise ``WebSocketClosed`` if the connection is closed."""
        if self.closed:
            raise WebSocketClosed("The connection is closed.")

    def _close_socket(self) -> None:
        """Close the socket and ignore errors."""
        try:
            self._sock.close()
        except OSError:
            pass

    def _abort(self) -> None:
        """Close the socket without a close frame: the stream is broken."""
        self.closed = True
        try:
            self._sock.close()
        except OSError:
            pass

    def close(self) -> None:
        """Close the connection. It is safe to call it more than one time."""
        if not self.closed:
            self.closed = True
            try:
                self._send_frame(OP_CLOSE, struct.pack(">H", 1000))
            except WebSocketError:
                pass
        try:
            self._sock.close()
        except OSError:
            pass
