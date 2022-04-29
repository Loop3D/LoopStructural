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
    logger.warning(
        "Cannot use MapView or RotationAnglePlotter as matplotlib is not installed. \n"
        "Install matplotlib and try again. "
    )
try:
    from .sphinx_scraper import _get_loop_visualisation_scraper
except:
    logger.error("Cannot use sphinx scraper, pip install -r docs/requirements.txt")
try:
    from .lavavu import LavaVuModelViewer
except ImportError:
    logger.error("Missing lavavu, can't import LavaVuModelViewer")

try:
    from .vtk_exporter import VtkExporter
except ImportError as e:
    logger.warning("Vtk export disabled: pip install meshio")
from ._scalar_field import ScalarField
from ._dash_view import DashView
