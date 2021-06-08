"""
Visualisation
=============

"""
from ..utils import getLogger
logger = getLogger(__name__)
try:
    import matplotlib
    from .map_viewer import MapView
    from .rotation_angle_plotter import RotationAnglePlotter
except ImportError:
    logger.warning('Cannot use MapView or RotationAnglePlotter as matplotlib is not installed. \n'\
                   'Install matplotlib and try again. ')
try:
    from .sphinx_scraper import _get_loop_visualisation_scraper
except:
    logger.error('Cannot use sphinx scraper, pip install -r docs/requirements.txt')
try:
    from .model_visualisation import LavaVuModelViewer
except:
    logger.error("Missing lavavu, can't import LavaVuModelViewer")

from ._scalar_field import ScalarField