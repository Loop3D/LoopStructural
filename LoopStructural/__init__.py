"""
LoopStructural
==============

"""

import logging
from logging.config import dictConfig

__all__ = ["GeologicalModel"]
from pathlib import Path
import tempfile

from .version import __version__

experimental = False
ch = logging.StreamHandler()
formatter = logging.Formatter("%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s")
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
loggers = {}
from .datatypes import BoundingBox
from .interpolators import InterpolatorBuilder
from .interpolators._api import LoopInterpolator
from .modelling.core.geological_model import GeologicalModel
from .utils import get_levels, getLogger, log_to_console, log_to_file, rng

logger = getLogger(__name__)
logger.info("Imported LoopStructural")


def setLogging(level="info"):
    """
    Set the logging parameters for log file

    Parameters
    ----------
    filename : string
        name of file or path to file
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    import LoopStructural

    logger = getLogger(__name__)

    levels = get_levels()
    level = levels.get(level, logging.WARNING)
    LoopStructural.ch.setLevel(level)

    for name in LoopStructural.loggers:
        logger = logging.getLogger(name)
        logger.setLevel(level)
    logger.info(f'Set logging to {level}')
