"""Tests for viewer discovery in ``loop_link.discovery``."""

from __future__ import annotations

import dataclasses
import json
import socket
import sys
from pathlib import Path

import pytest

from loop_link import discovery
from loop_link.discovery import DEFAULT_PORT, ViewerAddress, find_viewer, settings_dir
from loop_link.errors import NoViewerError, TokenError


@pytest.fixture(autouse=True)
def _no_env_token(monkeypatch):
    # The real environment must not change the result of a test.
    monkeypatch.delenv("GEOVIEWER_TOKEN", raising=False)


@pytest.fixture
def listeners():
    # Open sockets stay open until the end of the test.
    opened: list[socket.socket] = []

    def make() -> int:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind(("127.0.0.1", 0))
        sock.listen()
        opened.append(sock)
        return sock.getsockname()[1]

    yield make
    for sock in opened:
        sock.close()


def free_port() -> int:
    # Bind to port 0, read the port, close the socket: no program listens there.
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port


def write_link(settings: Path, token: str = "file-token") -> None:
    settings.mkdir(parents=True, exist_ok=True)
    (settings / "link.json").write_text(json.dumps({"token": token}), encoding="utf-8")


def write_instance(
    settings: Path, pid: int, port: int, kind: str = "window", focused: int = 0
) -> None:
    folder = settings / "instances"
    folder.mkdir(parents=True, exist_ok=True)
    data = {
        "version": 1,
        "pid": pid,
        "port": port,
        "kind": kind,
        "started": 1000,
        "focused": focused,
    }
    (folder / f"{pid}.json").write_text(json.dumps(data), encoding="utf-8")


def test_default_port_value():
    assert DEFAULT_PORT == 47615


def test_viewer_address_is_frozen():
    address = ViewerAddress(host="127.0.0.1", port=1, token="t")
    assert (address.host, address.port, address.token) == ("127.0.0.1", 1, "t")
    with pytest.raises(dataclasses.FrozenInstanceError):
        address.port = 2  # type: ignore[misc]


# Token order


def test_token_argument_wins(tmp_path, monkeypatch):
    monkeypatch.setenv("GEOVIEWER_TOKEN", "env-token")
    write_link(tmp_path, "file-token")
    address = find_viewer(port=1234, token="arg-token", settings=tmp_path)
    assert address.token == "arg-token"


def test_env_token_before_file(tmp_path, monkeypatch):
    monkeypatch.setenv("GEOVIEWER_TOKEN", "env-token")
    write_link(tmp_path, "file-token")
    address = find_viewer(port=1234, settings=tmp_path)
    assert address.token == "env-token"


def test_empty_env_token_is_not_a_token(tmp_path, monkeypatch):
    monkeypatch.setenv("GEOVIEWER_TOKEN", "")
    write_link(tmp_path, "file-token")
    address = find_viewer(port=1234, settings=tmp_path)
    assert address.token == "file-token"


def test_file_token_is_used(tmp_path):
    write_link(tmp_path, "file-token")
    address = find_viewer(port=1234, settings=tmp_path)
    assert address.token == "file-token"


def test_no_token_raises_token_error(tmp_path):
    with pytest.raises(TokenError):
        find_viewer(port=1234, settings=tmp_path)


def test_empty_env_and_no_file_raises_token_error(tmp_path, monkeypatch):
    monkeypatch.setenv("GEOVIEWER_TOKEN", "")
    with pytest.raises(TokenError):
        find_viewer(port=1234, settings=tmp_path)


def test_damaged_link_file_raises_token_error(tmp_path):
    (tmp_path / "link.json").write_text("this is not json", encoding="utf-8")
    with pytest.raises(TokenError):
        find_viewer(port=1234, settings=tmp_path)


def test_link_file_without_token_field_raises_token_error(tmp_path):
    (tmp_path / "link.json").write_text(json.dumps({"other": "value"}), encoding="utf-8")
    with pytest.raises(TokenError):
        find_viewer(port=1234, settings=tmp_path)


def test_link_file_with_empty_token_raises_token_error(tmp_path):
    write_link(tmp_path, "")
    with pytest.raises(TokenError):
        find_viewer(port=1234, settings=tmp_path)


# Port and host


def test_explicit_port_and_default_host(tmp_path):
    write_link(tmp_path)
    address = find_viewer(port=4321, settings=tmp_path)
    assert address == ViewerAddress(host="127.0.0.1", port=4321, token="file-token")


def test_explicit_host_is_kept(tmp_path):
    write_link(tmp_path)
    address = find_viewer(port=4321, host="localhost", settings=tmp_path)
    assert address.host == "localhost"


def test_explicit_loopback_host_is_returned(tmp_path):
    write_link(tmp_path)
    address = find_viewer(port=4321, host="127.0.0.1", settings=tmp_path)
    assert address.host == "127.0.0.1"


def test_explicit_port_skips_instance_files(tmp_path, listeners):
    write_link(tmp_path)
    live = listeners()
    write_instance(tmp_path, 111, live)
    # The explicit port has no listener. The function must not check it
    # and must not read the instance files.
    explicit = free_port()
    address = find_viewer(port=explicit, settings=tmp_path)
    assert address.port == explicit


# Instance files


def test_live_instance_file_is_found(tmp_path, listeners):
    write_link(tmp_path)
    live = listeners()
    write_instance(tmp_path, 111, live)
    address = find_viewer(settings=tmp_path)
    assert address.port == live
    assert address.host == "127.0.0.1"
    assert address.token == "file-token"


def test_file_without_listener_is_dropped(tmp_path, listeners):
    write_link(tmp_path)
    live = listeners()
    dead = free_port()
    # The dead file is better by kind and by focus time. It must still lose.
    write_instance(tmp_path, 111, dead, kind="window", focused=9000)
    write_instance(tmp_path, 222, live, kind="headless", focused=1)
    address = find_viewer(settings=tmp_path)
    assert address.port == live


def test_window_before_headless(tmp_path, listeners):
    write_link(tmp_path)
    window = listeners()
    headless = listeners()
    write_instance(tmp_path, 111, headless, kind="headless", focused=9000)
    write_instance(tmp_path, 222, window, kind="window", focused=1)
    address = find_viewer(settings=tmp_path)
    assert address.port == window


def test_largest_focused_wins_in_same_kind(tmp_path, listeners):
    write_link(tmp_path)
    old = listeners()
    recent = listeners()
    write_instance(tmp_path, 111, recent, kind="window", focused=500)
    write_instance(tmp_path, 222, old, kind="window", focused=100)
    address = find_viewer(settings=tmp_path)
    assert address.port == recent


def test_largest_focused_wins_for_headless(tmp_path, listeners):
    write_link(tmp_path)
    old = listeners()
    recent = listeners()
    write_instance(tmp_path, 111, old, kind="headless", focused=100)
    write_instance(tmp_path, 222, recent, kind="headless", focused=500)
    address = find_viewer(settings=tmp_path)
    assert address.port == recent


def test_damaged_instance_file_is_ignored(tmp_path, listeners):
    write_link(tmp_path)
    live = listeners()
    write_instance(tmp_path, 222, live)
    (tmp_path / "instances" / "111.json").write_text("not json", encoding="utf-8")
    (tmp_path / "instances" / "333.json").write_text(json.dumps({"pid": 333}), encoding="utf-8")
    address = find_viewer(settings=tmp_path)
    assert address.port == live


# No viewer


def test_no_file_and_no_listener_raises_no_viewer(tmp_path, monkeypatch):
    write_link(tmp_path)
    # Use a free port as default port: the real machine does not change the result.
    monkeypatch.setattr(discovery, "DEFAULT_PORT", free_port())
    with pytest.raises(NoViewerError):
        find_viewer(settings=tmp_path)


def test_only_dead_files_raises_no_viewer(tmp_path, monkeypatch):
    write_link(tmp_path)
    monkeypatch.setattr(discovery, "DEFAULT_PORT", free_port())
    write_instance(tmp_path, 111, free_port())
    with pytest.raises(NoViewerError):
        find_viewer(settings=tmp_path)


def test_default_port_listener_is_used_without_files(tmp_path, monkeypatch, listeners):
    write_link(tmp_path)
    live = listeners()
    monkeypatch.setattr(discovery, "DEFAULT_PORT", live)
    address = find_viewer(settings=tmp_path)
    assert address.port == live


# settings_dir


def test_settings_dir_windows(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setenv("APPDATA", str(tmp_path))
    assert settings_dir() == tmp_path / "geoviewer"


def test_settings_dir_other_platform(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setenv("HOME", str(tmp_path))
    # Path.home() reads USERPROFILE on Windows.
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    assert settings_dir() == tmp_path / ".geoviewer"


def test_non_finite_focused_counts_as_zero(tmp_path, listeners):
    write_link(tmp_path)
    nan_port = listeners()
    good_port = listeners()
    write_instance(tmp_path, 111, good_port, kind="window", focused=5)
    text = json.dumps(
        {"version": 1, "pid": 222, "port": nan_port, "kind": "window", "started": 1000}
    )
    text = text[:-1] + ', "focused": NaN}'
    (tmp_path / "instances" / "222.json").write_text(text, encoding="utf-8")
    address = find_viewer(settings=tmp_path)
    assert address.port == good_port
