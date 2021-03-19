"""


LoopStructural API
=======================

.. autosummary::
    :toctree:


    interpolators
    modelling
    visualisation
    utils

"""

import logging
from logging.config import dictConfig
import tempfile
from pathlib import Path
#set up logging
# temp_file = tempfile.mkdtemp()
# if temp_file:
#     # temp_file = tempfile.tempdir+Path('/default-loop-structural-logfile.log')
#     log_to_file(temp_file)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s ~ %(name)-12s ~ %(levelname)-10s ~ %(message)s')
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
loggers = {}
__version__ = '1.0.9'
from .modelling.core.geological_model import GeologicalModel
from .utils import log_to_console, log_to_file, getLogger
logger = getLogger(__name__)
logger.info("Imported LoopStructural")
