"""
LoopStructural
==============

"""

import logging

__all__ = ["GeologicalModel"]

experimental = False
ch = logging.StreamHandler()
formatter = logging.Formatter("%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s")
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
loggers = {}
from .modelling.core.geological_model import GeologicalModel
from .utils import getLogger

logger = getLogger(__name__)
logger.info("Imported LoopStructural")
