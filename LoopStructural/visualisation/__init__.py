from ..utils import getLogger

logger = getLogger(__name__)

try:
    from loopstructuralvisualisation import (
        Loop3DView,
        RotationAnglePlotter,
        Loop2DView,
        StratigraphicColumnView,
    )
except ImportError as e:
    logger.error("Please install the loopstructuralvisualisation package")
    logger.error("pip install loopstructuralvisualisation")
    raise e
