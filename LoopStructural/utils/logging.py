import logging
import os
from typing import Dict, Optional, Union

import LoopStructural
from loop_common.logging import (
    LogSink,
    StreamSink,
    FileSink,
    SqliteSink,
    timed_stage,
    timed,
)
from loop_common.logging.sinks import LogCallable, _CallableHandler

from ._api_registry import public_api

__all__ = [
    "getLogger",
    "log_to_file",
    "log_to_console",
    "get_levels",
    "LogSink",
    "StreamSink",
    "FileSink",
    "SqliteSink",
    "add_sink",
    "remove_sink",
    "timed_stage",
    "timed",
]


def get_levels():
    """dict for converting to logger levels from string


    Returns
    -------
    dict
        contains all strings with corresponding logging levels.
    """
    return {
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "debug": logging.DEBUG,
    }


@public_api(tier="stable")
def getLogger(name):
    """Get (or create) a stdlib `logging.Logger` wired into LoopStructural's shared sinks.

    The returned object is a genuine `logging.Logger`, so host applications
    (e.g. the QGIS plugin) can keep attaching their own handlers to it
    directly, exactly as before. `LoopStructural.utils.add_sink` is the
    higher-level, documented way to do the same thing -- as a `LogSink`
    subclass or a plain callable -- without reaching into stdlib logging
    internals, and without needing to re-attach to loggers created later.

    Parameters
    ----------
    name : str
        Logger name, conventionally `__name__` of the calling module.

    Returns
    -------
    logging.Logger
    """
    logger = logging.getLogger(name)
    logger.addHandler(LoopStructural.ch)
    for handler in LoopStructural._extra_sinks:
        logger.addHandler(handler)
    # don't pass message back up the chain, what an odd default behavior
    logger.propagate = False
    # store the loopstructural loggers so we can change values
    LoopStructural.loggers[name] = logger
    return logger


def log_to_file(filename, overwrite=True, level="info"):
    """Set the logging parameters for log file


    Parameters
    ----------
    filename : string
        name of file or path to file
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    logger = getLogger(__name__)
    if os.path.isfile(filename):
        logger.warning("Overwriting existing logfile. To avoid this, set overwrite=False")
        os.remove(filename)
    levels = get_levels()
    level = levels.get(level, logging.WARNING)
    fh = logging.FileHandler(filename)
    fh.setFormatter(LoopStructural.formatter)
    fh.setLevel(level)
    for logger in LoopStructural.loggers.values():
        for hdlr in logger.handlers[:]:  # remove the existing file handlers
            if isinstance(hdlr, logging.FileHandler):  # fixed two typos here
                logger.removeHandler(hdlr)
        logger.addHandler(fh)
        logger.setLevel(level)


def log_to_console(level="warning"):
    """Set the level of logging to the console


    Parameters
    ----------
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level, logging.WARNING)
    for logger in LoopStructural.loggers.values():
        for hdlr in logger.handlers:
            # both stream and file are base stream, so check if not a filehandler
            if not isinstance(hdlr, logging.FileHandler):
                logger.removeHandler(hdlr)
                hdlr = LoopStructural.ch
                hdlr.setLevel(level)
                logger.addHandler(hdlr)


@public_api(tier="provisional")
def add_sink(
    sink: Union[LogSink, LogCallable], *, loggers: Optional[Dict[str, logging.Logger]] = None
) -> logging.Handler:
    """Attach a sink to every currently-registered LoopStructural logger.

    Parameters
    ----------
    sink : LogSink | Callable[[logging.LogRecord], None]
        A `LogSink` subclass instance, or a plain callable -- both are
        supported extension points for host applications (see `LogSink`).
    loggers : dict[str, logging.Logger], optional
        Registry to attach to; defaults to `LoopStructural.loggers`.

    Returns
    -------
    logging.Handler
        The resulting handler, so it can later be detached with `remove_sink`.

    Notes
    -----
    Loggers created with `getLogger` *after* this call also pick up the
    sink automatically, matching how the built-in console sink already
    behaves.
    """
    handler = sink.handler() if isinstance(sink, LogSink) else _CallableHandler(sink)
    LoopStructural._extra_sinks.append(handler)
    target = loggers if loggers is not None else LoopStructural.loggers
    for logger in target.values():
        logger.addHandler(handler)
    return handler


@public_api(tier="provisional")
def remove_sink(
    handler: logging.Handler, *, loggers: Optional[Dict[str, logging.Logger]] = None
) -> None:
    """Detach a handler previously returned by `add_sink`."""
    if handler in LoopStructural._extra_sinks:
        LoopStructural._extra_sinks.remove(handler)
    target = loggers if loggers is not None else LoopStructural.loggers
    for logger in target.values():
        logger.removeHandler(handler)
