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
from .modelling.core.geological_model import GeologicalModel
from .visualisation.model_visualisation import LavaVuModelViewer
from .visualisation.map_viewer import MapView
from .utils.utils import log_to_console, log_to_file

log_to_file('/tmp/default-loop-structural-logfile.log')
log_to_console()