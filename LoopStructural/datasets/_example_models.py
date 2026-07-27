from ..utils import getLogger

logger = getLogger(__name__)

vis = True
try:
    pass
except Exception:
    logger.warning("No visualisation")
    vis = False


def _build_claudius():
    pass
