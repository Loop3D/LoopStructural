"""Find a running viewer and the token to talk to it.

The functions here only read settings files and test a socket. They never
start a program.
"""

from __future__ import annotations

import ipaddress
import json
import math
import os
import socket
import sys
import warnings
from dataclasses import dataclass, field
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
    token: str = field(repr=False)


def _is_loopback(host: str) -> bool:
    if host.lower() == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


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


def _resolve_token(token: str | None, settings: Path, local: bool) -> str:
    if token:
        return token
    env_token = os.environ.get("GEOVIEWER_TOKEN")
    if env_token:
        return env_token
    # The saved token belongs to the viewer on this computer: never send it to another host.
    file_token = _file_token(settings) if local else ""
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
        # The instance files list the ports of the viewers on this computer only.
        files = sorted(folder.glob("*.json")) if _is_loopback(host) else []
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
        ``GEOVIEWER_TOKEN`` environment variable, then ``link.json``. For a
        host that is not loopback, it does not read ``link.json``.
    host : str
        The host of the viewer. For a host that is not loopback, the function
        does not read the instance files and gives a warning, because the link
        is not encrypted.
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
    local = _is_loopback(host)
    if not local:
        warnings.warn(
            f"The link to {host} is not encrypted: the token and the data are sent in plain text.",
            UserWarning,
            stacklevel=2,
        )
    resolved = _resolve_token(token, folder, local)
    if port is not None:
        return ViewerAddress(host=host, port=port, token=resolved)
    return ViewerAddress(host=host, port=_find_port(folder, host, timeout), token=resolved)
