"""
Analysis
========

Various tools for analysing loopstructural models, including calculating fault intersections and fault toplogies
"""
from LoopStructural.utils import getLogger
import LoopStructural

logger = getLogger(__name__)
if LoopStructural.experimental:
    logger.warning(
        "LoopStructural.analysis is experimental and may not perform as expected"
    )
    from ._fault_displacement import displacement_missfit
    from ._fault_intersection import calculate_fault_intersections
    from ._topology import calculate_fault_topology_matrix
