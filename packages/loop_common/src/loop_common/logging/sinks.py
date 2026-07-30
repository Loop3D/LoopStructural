"""Pluggable log-sink infrastructure, shared across Loop packages.

:class:`LogSink` is the documented extension point host applications use to
route a package's log records into their own systems -- either by
subclassing it, or by passing a plain
``Callable[[logging.LogRecord], None]`` straight to a package's
``add_sink`` (e.g. ``LoopStructural.utils.add_sink``), no subclassing
required.

Attaching sinks to a specific logger registry (LoopStructural keeps its own
in ``LoopStructural.loggers``/``LoopStructural._extra_sinks``) is the
caller's responsibility -- these classes only build the
``logging.Handler`` that gets attached.
"""

from __future__ import annotations

import logging
import sqlite3
import threading
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Union

LogCallable = Callable[[logging.LogRecord], None]

__all__ = [
    "LogSink",
    "StreamSink",
    "FileSink",
    "SqliteSink",
    "default_formatter",
]


def default_formatter() -> logging.Formatter:
    """Return the formatter used by the built-in sinks."""
    return logging.Formatter("%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s")


class LogSink(ABC):
    """Base class for a pluggable logging destination.

    Subclass and implement :meth:`emit` to receive every
    ``logging.LogRecord`` forwarded to a logger. A plain
    ``Callable[[logging.LogRecord], None]`` works too and does not require
    subclassing this class at all.
    """

    level: int = logging.NOTSET

    @abstractmethod
    def emit(self, record: logging.LogRecord) -> None:
        """Handle a single log record."""

    def handler(self) -> logging.Handler:
        """Build the ``logging.Handler`` used to attach this sink to a logger."""
        return _CallableHandler(self.emit, level=self.level)


class _CallableHandler(logging.Handler):
    """Adapts a plain callable (or a `LogSink.emit`) to `logging.Handler`."""

    def __init__(self, callback: LogCallable, *, level: int = logging.NOTSET):
        super().__init__(level=level)
        self._callback = callback

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._callback(record)
        except Exception:
            self.handleError(record)


class StreamSink(LogSink):
    """Writes formatted records to a stream, defaulting to stderr."""

    def __init__(
        self,
        stream=None,
        *,
        formatter: Optional[logging.Formatter] = None,
        level: int = logging.WARNING,
    ):
        self.level = level
        self._handler = logging.StreamHandler(stream)
        self._handler.setFormatter(formatter or default_formatter())
        self._handler.setLevel(level)

    def emit(self, record: logging.LogRecord) -> None:
        self._handler.emit(record)

    def handler(self) -> logging.Handler:
        return self._handler


class FileSink(LogSink):
    """Writes formatted records to a log file, creating parent directories as needed."""

    def __init__(
        self,
        path: Union[str, Path],
        *,
        overwrite: bool = False,
        formatter: Optional[logging.Formatter] = None,
        level: int = logging.INFO,
    ):
        self.path = Path(path)
        self.level = level
        if overwrite and self.path.exists():
            self.path.unlink()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handler = logging.FileHandler(self.path)
        self._handler.setFormatter(formatter or default_formatter())
        self._handler.setLevel(level)

    def emit(self, record: logging.LogRecord) -> None:
        self._handler.emit(record)

    def handler(self) -> logging.Handler:
        return self._handler


class SqliteSink(LogSink):
    """Writes structured log records to a SQLite database for querying run history.

    Records produced by ``timed_stage``/``timed`` carry extra attributes
    (``stage``, ``event``, ``duration_s``, ``run_id``) which are stored in
    dedicated columns, so build/interpolation timings can be queried
    directly (``sink.query(stage="update")``) instead of parsed out of
    formatted log text.
    """

    _COLUMNS = (
        "timestamp",
        "logger_name",
        "level",
        "message",
        "module",
        "func_name",
        "lineno",
        "stage",
        "event",
        "duration_s",
        "run_id",
    )

    def __init__(
        self, path: Union[str, Path], *, table: str = "log_records", level: int = logging.NOTSET
    ):
        self.path = Path(path)
        self.level = level
        self.table = table
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        with self._lock, self._connect() as conn:
            columns_sql = ", ".join(
                f"{c} REAL" if c == "duration_s" else f"{c} TEXT" for c in self._COLUMNS
            )
            conn.execute(
                f"CREATE TABLE IF NOT EXISTS {self.table} "
                f"(id INTEGER PRIMARY KEY AUTOINCREMENT, {columns_sql})"
            )

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, check_same_thread=False)

    def emit(self, record: logging.LogRecord) -> None:
        row = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "logger_name": record.name,
            "level": record.levelname,
            "message": record.getMessage(),
            "module": record.module,
            "func_name": record.funcName,
            "lineno": record.lineno,
            "stage": getattr(record, "stage", None),
            "event": getattr(record, "event", None),
            "duration_s": getattr(record, "duration_s", None),
            "run_id": getattr(record, "run_id", None),
        }
        columns = ", ".join(row)
        placeholders = ", ".join("?" for _ in row)
        with self._lock, self._connect() as conn:
            conn.execute(
                f"INSERT INTO {self.table} ({columns}) VALUES ({placeholders})",
                tuple(row.values()),
            )

    def query(
        self,
        *,
        stage: Optional[str] = None,
        run_id: Optional[str] = None,
        logger_name: Optional[str] = None,
        level: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> List[Dict]:
        """Query recorded log rows, optionally filtered. Returns dict rows, oldest first."""
        clauses, params = [], []
        for column, value in (
            ("stage", stage),
            ("run_id", run_id),
            ("logger_name", logger_name),
            ("level", level),
        ):
            if value is not None:
                clauses.append(f"{column} = ?")
                params.append(value)
        sql = f"SELECT * FROM {self.table}"
        if clauses:
            sql += " WHERE " + " AND ".join(clauses)
        sql += " ORDER BY id"
        if limit is not None:
            sql += " LIMIT ?"
            params.append(int(limit))
        with self._lock, self._connect() as conn:
            conn.row_factory = sqlite3.Row
            rows = conn.execute(sql, params).fetchall()
        return [dict(row) for row in rows]
