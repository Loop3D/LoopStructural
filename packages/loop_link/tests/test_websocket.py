"""Tests for the small WebSocket client in ``loop_link._websocket``."""

from __future__ import annotations

import socket
import struct
import threading
import time

import pytest

from loop_link._websocket import OP_BINARY, WebSocketClient


def test_large_frame_is_read_in_linear_time() -> None:
    # Reads of a growing bytes buffer were quadratic: 16 MB took seconds.
    size = 32 * 1024 * 1024
    ours, theirs = socket.socketpair()
    head = bytes([0x80 | OP_BINARY, 127]) + struct.pack(">Q", size)
    writer = threading.Thread(target=theirs.sendall, args=(head + bytes(size),))
    writer.start()
    try:
        start = time.monotonic()
        data = WebSocketClient(ours).receive(10.0)
        elapsed = time.monotonic() - start
    finally:
        writer.join()
        ours.close()
        theirs.close()
    assert data is not None and len(data) == size
    assert elapsed < 2.0


class _InterruptedSocket:
    """A socket whose second sendall is interrupted, after part of the frame went out."""

    def __init__(self) -> None:
        self.calls = 0
        self.closed = False

    def settimeout(self, value: float | None) -> None:
        pass

    def sendall(self, data: bytes) -> None:
        self.calls += 1
        if self.calls == 2:
            raise KeyboardInterrupt

    def close(self) -> None:
        self.closed = True


def test_interrupt_during_send_closes_the_client() -> None:
    sock = _InterruptedSocket()
    client = WebSocketClient(sock)  # type: ignore[arg-type]
    with pytest.raises(KeyboardInterrupt):
        client.send_binary(b"x" * 1000)
    assert client.closed and sock.closed


def test_interrupt_while_waiting_for_a_frame_keeps_the_client() -> None:
    # No byte of the frame was read, so the stream is still aligned.
    ours, theirs = socket.socketpair()
    real = ours

    class InterruptedOnce:
        interrupted = False

        def __getattr__(self, name: str) -> object:
            return getattr(real, name)

        def recv(self, size: int) -> bytes:
            if not InterruptedOnce.interrupted:
                InterruptedOnce.interrupted = True
                raise KeyboardInterrupt
            return real.recv(size)

    client = WebSocketClient(InterruptedOnce())  # type: ignore[arg-type]
    try:
        with pytest.raises(KeyboardInterrupt):
            client.receive(1.0)
        assert not client.closed
        theirs.sendall(bytes([0x80 | OP_BINARY, 5]) + b"hello")
        assert client.receive(1.0) == b"hello"
    finally:
        ours.close()
        theirs.close()
