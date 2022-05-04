"""
Geological modelling classes and functions

"""
__all__ = ['GeologicalModel', 'ProcessInputData', 'LavaVuModelViewer', 'Map2LoopProcessor','LoopProjectfileProcessor']
from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError

logger = getLogger(__name__)
from LoopStructural.modelling.input import (
    ProcessInputData,
    Map2LoopProcessor,
)

try:
    from LoopStructural.modelling.input.project_file import LoopProjectfileProcessor
except (LoopImportError, ImportError):
    logger.warning(
        "Cannot use LoopProjectfileProcessor: Loop project file cannot be imported, try installing LoopProjectFile"
    )
# from LoopStructural.modelling.features import StructuralFrame
# from LoopStructural.modelling.features.fault import FaultSegment
