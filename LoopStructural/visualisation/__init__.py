from ..utils import getLogger

logger = getLogger(__name__)

try:
    from loopstructuralvisualisation import (
        Loop2DView,
        Loop3DView,
        RotationAnglePlotter,
        StratigraphicColumnView,
    )
except ImportError:
    logger.error("Please install the loopstructuralvisualisation package")
    logger.error("pip install loopstructuralvisualisation")
    raise
