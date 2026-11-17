"""Find a running viewer and the token to talk to it.

The functions here only read settings files and test a socket. They never
start a program.
"""

from __future__ import annotations

import json
import math
import os
import socket
import sys
from dataclasses import dataclass
from pathlib import Path

from loop_link.errors import NoViewerError, TokenError

DEFAULT_PORT: int = 47615

_TOKEN_TEXT = "Set GEOVIEWER_TOKEN, or save a token in the link popup of the viewer."
_NO_VIEWER_TEXT = (
    "No viewer found. Start geoviewer and turn the link on. " "This library never starts a viewer."
)


@dataclass(frozen=True)
class ViewerAddress:
    """The address and token of a viewer.

    Parameters
    ----------
    host : str
        The host name or IP address of the viewer.
    port : int
        The TCP port of the viewer.
    token : str
        The access token of the viewer.
    """

    host: str
    port: int
    token: str


def settings_dir() -> Path:
    r"""Return the folder where the viewer keeps its settings.

    Returns
    -------
    pathlib.Path
        ``%APPDATA%\\geoviewer`` on Windows, else ``~/.geoviewer``.
    """
    if sys.platform == "win32":
        appdata = os.environ.get("APPDATA")
        base = Path(appdata) if appdata else Path.home() / "AppData" / "Roaming"
        return base / "geoviewer"
    return Path.home() / ".geoviewer"


def _read_json(path: Path) -> object:
    # Return the parsed JSON, or None if the file is missing or damaged.
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _file_token(settings: Path) -> str:
    data = _read_json(settings / "link.json")
    if isinstance(data, dict):
        token = data.get("token")
        if isinstance(token, str):
            return token
    return ""


def _resolve_token(token: str | None, settings: Path) -> str:
    if token:
        return token
    env_token = os.environ.get("GEOVIEWER_TOKEN")
    if env_token:
        return env_token
    file_token = _file_token(settings)
    if file_token:
        return file_token
    raise TokenError(_TOKEN_TEXT)


def _is_alive(host: str, port: int, timeout: float) -> bool:
    try:
        sock = socket.create_connection((host, port), timeout)
    except OSError:
        return False
    sock.close()
    return True


def _find_port(settings: Path, host: str, timeout: float) -> int:
    live: list[tuple[tuple[bool, int], int]] = []
    folder = settings / "instances"
    try:
        files = sorted(folder.glob("*.json"))
    except OSError:
        files = []
    for path in files:
        data = _read_json(path)
        if not isinstance(data, dict):
            continue
        port = data.get("port")
        kind = data.get("kind")
        focused = data.get("focused", 0)
        if isinstance(focused, bool) or not isinstance(focused, (int, float)):
            focused = 0
        if not math.isfinite(focused):
            focused = 0
        if not isinstance(port, int) or isinstance(port, bool) or not isinstance(kind, str):
            continue
        if _is_alive(host, port, timeout):
            live.append(((kind != "window", -focused), port))
    if live:
        live.sort(key=lambda item: item[0])
        return live[0][1]
    if _is_alive(host, DEFAULT_PORT, timeout):
        return DEFAULT_PORT
    raise NoViewerError(_NO_VIEWER_TEXT)


def find_viewer(
    *,
    port: int | None = None,
    token: str | None = None,
    host: str = "127.0.0.1",
    settings: Path | None = None,
    timeout: float = 0.5,
) -> ViewerAddress:
    """Find a running viewer and the token for it.

    Parameters
    ----------
    port : int, optional
        The port of the viewer. If given, it is returned as given: the
        function does not check it and does not read the instance files.
    token : str, optional
        The access token. If not given, the function reads the
        ``GEOVIEWER_TOKEN`` environment variable, then ``link.json``.
    host : str
        The host of the viewer.
    settings : pathlib.Path, optional
        The settings folder. The default is :func:`settings_dir`.
    timeout : float
        The time in seconds to wait for each connection test.

    Returns
    -------
    ViewerAddress
        The host, port and token of the viewer.

    Raises
    ------
    TokenError
        If no token is available.
    NoViewerError
        If no viewer is running.
    """
    folder = settings if settings is not None else settings_dir()
    resolved = _resolve_token(token, folder)
    if port is not None:
        return ViewerAddress(host=host, port=port, token=resolved)
    return ViewerAddress(host=host, port=_find_port(folder, host, timeout), token=resolved)
