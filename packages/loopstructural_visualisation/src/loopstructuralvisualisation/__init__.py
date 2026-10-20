from LoopStructural.utils import getLogger

from ._3d_viewer import Loop3DView
from ._rotation_angle import RotationAnglePlotter
from ._2d_viewer import Loop2DView
from ._stratigraphic_column import StratigraphicColumnView

logger = getLogger(__name__)

try:
    from . import _register_loop_ui  # noqa: F401  # registers the optional notebook UI
except ImportError:
    logger.warning("Could not import trame ui")
