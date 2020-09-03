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
from .modelling.core.geological_model import GeologicalModel
from .visualisation.model_visualisation import LavaVuModelViewer
from .visualisation.map_viewer import MapView
from .utils.utils import log_to_console, log_to_file

log_to_console()
__version__ = '0.0.1'
