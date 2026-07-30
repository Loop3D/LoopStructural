"""Timing/instrumentation helpers for staged, long-running work.

:func:`timed_stage` is the primitive (a context manager); :func:`timed` is
a thin decorator wrapping it for whole-function timing. Both emit
structured start/end log records (via the `extra=` mechanism of the
stdlib `logging` module) carrying `stage`, `event`, `run_id` and, on
completion, `duration_s` -- fields a `SqliteSink` stores in dedicated
columns so run history can be queried without parsing message text.
"""

from __future__ import annotations

import functools
import logging
import time
import uuid
from contextlib import contextmanager
from typing import Callable, Optional

__all__ = ["timed_stage", "timed"]


@contextmanager
def timed_stage(
    logger: logging.Logger,
    stage: str,
    *,
    run_id: Optional[str] = None,
    level: int = logging.INFO,
    **extra,
):
    """Time a named stage (e.g. "update", "interpolate") and log its duration.

    Emits a ``event="start"`` record on entry and an ``event="end"`` record
    (with a ``duration_s`` field) on exit -- even if the block raises.

    Parameters
    ----------
    logger : logging.Logger
        Logger to emit the start/end records on.
    stage : str
        Name of the stage being timed, e.g. "update" or "interpolate".
    run_id : str, optional
        Correlates stages from the same run; generated if omitted.
    level : int, optional
        Logging level for the emitted records, by default `logging.INFO`.
    **extra
        Additional fields attached to both log records.

    Yields
    ------
    str
        The `run_id` used for this timed block.
    """
    run_id = run_id or uuid.uuid4().hex[:8]
    start = time.perf_counter()
    # stacklevel=3: past this frame and contextlib's generator-CM __enter__,
    # so module/funcName/lineno on the record point at the `with` site.
    logger.log(
        level,
        f"{stage}: started",
        extra={"stage": stage, "event": "start", "run_id": run_id, **extra},
        stacklevel=3,
    )
    try:
        yield run_id
    finally:
        duration_s = time.perf_counter() - start
        logger.log(
            level,
            f"{stage}: finished in {duration_s:.3f}s",
            extra={
                "stage": stage,
                "event": "end",
                "run_id": run_id,
                "duration_s": duration_s,
                **extra,
            },
            stacklevel=3,
        )


def timed(
    stage: Optional[str] = None,
    *,
    logger: Optional[logging.Logger] = None,
    level: int = logging.INFO,
):
    """Decorator version of `timed_stage`, timing an entire function call.

    Parameters
    ----------
    stage : str, optional
        Name of the stage; defaults to the wrapped function's qualified name.
    logger : logging.Logger, optional
        Logger to use; defaults to a stdlib logger named after the
        function's module.
    level : int, optional
        Logging level for the emitted records, by default `logging.INFO`.
    """

    def decorator(func: Callable) -> Callable:
        stage_name = stage or func.__qualname__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            active_logger = logger or logging.getLogger(func.__module__)
            with timed_stage(active_logger, stage_name):
                return func(*args, **kwargs)

        return wrapper

    return decorator
