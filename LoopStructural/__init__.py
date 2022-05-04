"""
LoopStructural API
=======================

"""

import logging
from logging.config import dictConfig
__all__ = ['GeologicalModel']
import tempfile
from pathlib import Path
from .version import __version__

experimental = False
ch = logging.StreamHandler()
formatter = logging.Formatter(
    "%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s"
)
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
loggers = {}
from .modelling.core.geological_model import GeologicalModel
from .utils import log_to_console, log_to_file, getLogger

logger = getLogger(__name__)
logger.info("Imported LoopStructural")
