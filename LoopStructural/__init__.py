"""
LoopStructural
==============

"""

import logging
from dataclasses import dataclass
from logging.config import dictConfig

__all__ = [
    "BoundingBox",
    "FaultTopology",
    "FileSink",
    "GeologicalModel",
    "InterpolatorBuilder",
    "LogSink",
    "LoopInterpolator",
    "LoopStructuralConfig",
    "SqliteSink",
    "StratigraphicColumn",
    "StreamSink",
    "add_sink",
    "getLogger",
    "get_levels",
    "log_to_console",
    "log_to_file",
    "remove_sink",
    "rng",
    "setLogging",
    "timed",
    "timed_stage",
]
import tempfile
from pathlib import Path

from .version import __version__

experimental = False
ch = logging.StreamHandler()
formatter = logging.Formatter("%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s")
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
loggers = {}
# Handlers attached via LoopStructural.utils.add_sink(); replayed onto every
# logger getLogger() creates from here on, in addition to the default `ch`.
_extra_sinks = []


@dataclass
class LoopStructuralConfig:
    """Configuration for LoopStructural package.

    This dataclass holds configuration parameters for the LoopStructural
    geological modelling package.

    Parameters
    ----------
    nelements : int, optional
        The default number of elements to use in interpolation, by default 10_000

    Examples
    --------
    >>> config = LoopStructuralConfig(nelements=50000)
    >>> config.nelements
    50000
    """
   
    nelements: int = 10_000

from .geometry import BoundingBox
from .interpolators import InterpolatorBuilder
from .interpolators._api import LoopInterpolator
from .modelling.core.fault_topology import FaultTopology
from .modelling.core.geological_model import GeologicalModel
from .modelling.core.stratigraphic_column import StratigraphicColumn
from .utils import (
    FileSink,
    LogSink,
    SqliteSink,
    StreamSink,
    add_sink,
    get_levels,
    getLogger,
    log_to_console,
    log_to_file,
    remove_sink,
    rng,
    timed,
    timed_stage,
)

logger = getLogger(__name__)
logger.info("Imported LoopStructural")


def setLogging(level="info", handler=None):
    """Set the logging parameters for log file or custom handler.

    Parameters
    ----------
    level : str, optional
        Logging level to set, by default "info"
        Valid options: 'info', 'warning', 'error', 'debug'
    handler : logging.Handler, optional
        A logging handler to use instead of the default StreamHandler, by default None

    Examples
    --------
    >>> import LoopStructural
    >>> LoopStructural.setLogging('debug')
    >>> LoopStructural.setLogging('info', logging.FileHandler('loop.log'))
    """
    import LoopStructural

    levels = get_levels()
    level_value = levels.get(level, logging.WARNING)

    # Create default handler if none provided
    if handler is None:
        handler = logging.StreamHandler()

    formatter = logging.Formatter(
        "%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s"
    )
    handler.setFormatter(formatter)
    handler.setLevel(level_value)

    # Replace handlers in all known loggers
    for name in LoopStructural.loggers:
        logger = logging.getLogger(name)
        logger.handlers = []
        logger.addHandler(handler)
        logger.setLevel(level_value)

    # Also apply to main module logger
    main_logger = logging.getLogger(__name__)
    main_logger.handlers = []
    main_logger.addHandler(handler)
    main_logger.setLevel(level_value)

    main_logger.info(f"Set logging to {level}")
