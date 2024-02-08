"""
Analysis
========

Various tools for analysing loopstructural models, including calculating fault intersections and fault toplogies
"""
from ..utils import getLogger
import LoopStructural

logger = getLogger(__name__)
if LoopStructural.experimental:
    logger.warning("LoopStructural.analysis is experimental and may not perform as expected")
