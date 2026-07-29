"""Zero-boilerplate logging interface.

Usage::

    from lgutils.logging import get_logger

    log = get_logger(__name__)
    log.info("hello world")

    # With a log file:
    log = get_logger(__name__, log_file="run.log")
    log.warning("this also goes to run.log")
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Union

# ---------------------------------------------------------------------------
# Optional loguru detection
# ---------------------------------------------------------------------------

try:
    from loguru import logger as _loguru_logger  # noqa: F401

    _LOGURU_AVAILABLE = True
except ImportError:
    _LOGURU_AVAILABLE = False

_DEFAULT_FMT = "%(asctime)s | %(name)-20s | %(levelname)-8s | %(message)s"
_LOGURU_FMT = "{time:YYYY-MM-DD HH:mm:ss} | {name:<20} | {level:<8} | {message}"


def get_logger(
    name: str,
    level: Union[str, int] = "INFO",
    log_file: Union[str, Path, None] = None,
    fmt: Union[str, None] = None,
    use_loguru: Union[bool, None] = None,
):
    """Return a configured logger with no boilerplate required at the call site.

    Parameters
    ----------
    name:
        Logger name, typically ``__name__``.
    level:
        Log level string (``"DEBUG"``, ``"INFO"``, ``"WARNING"``, ``"ERROR"``)
        or the corresponding integer constant.
    log_file:
        Optional path to a log file.  Output is written to both stdout/stderr
        **and** the file.  Parent directories are created automatically.
    fmt:
        Custom format string.  For stdlib loggers this is a ``%``-style
        format; for loguru it is a loguru format string.
    use_loguru:
        Override auto-detection.  ``True`` forces loguru (raises
        ``RuntimeError`` if not installed).  ``False`` forces stdlib.
        ``None`` (default) uses loguru when available, stdlib otherwise.

    Returns
    -------
    A logger object with ``.debug``, ``.info``, ``.warning``,
    ``.error``, and ``.exception`` methods.
    """
    _use_loguru = _LOGURU_AVAILABLE if use_loguru is None else use_loguru
    if _use_loguru:
        return _build_loguru_logger(name, level, log_file, fmt)
    return _build_stdlib_logger(name, level, log_file, fmt)


# ---------------------------------------------------------------------------
# Loguru backend
# ---------------------------------------------------------------------------


def _build_loguru_logger(name, level, log_file, fmt):
    if not _LOGURU_AVAILABLE:
        raise RuntimeError("loguru is not installed. Install it with: pip install loguru")
    from loguru import logger

    # Remove the default handler so we configure our own sinks.
    logger.remove()
    _fmt = fmt or _LOGURU_FMT
    logger.add(sys.stderr, level=level, format=_fmt, colorize=True)
    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        logger.add(str(log_file), level=level, format=_fmt, enqueue=True)

    return logger.bind(name=name)


# ---------------------------------------------------------------------------
# Stdlib backend
# ---------------------------------------------------------------------------


def _build_stdlib_logger(name, level, log_file, fmt):
    log = logging.getLogger(name)
    log.setLevel(level)

    # Guard against duplicate handlers on repeated calls with the same name.
    if log.handlers:
        return log

    _fmt = fmt or _DEFAULT_FMT
    formatter = logging.Formatter(_fmt)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)

    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)
        log.addHandler(file_handler)

    log.propagate = False
    return log
