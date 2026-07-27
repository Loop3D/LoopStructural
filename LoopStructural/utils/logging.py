import logging
import LoopStructural
import os

from ._api_registry import public_api
from ._log_sinks import (
    LogSink,
    StreamSink,
    FileSink,
    SqliteSink,
    add_sink,
    remove_sink,
)
from ._log_timing import timed_stage, timed

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
