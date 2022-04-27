"""
Geological modelling classes and functions

"""
from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError

logger = getLogger(__name__)
from LoopStructural.modelling.input.process_data import ProcessInputData

try:
    from LoopStructural.modelling.input.project_file import LoopProjectfileProcessor
except (LoopImportError, ImportError):
    logger.warning(
        "Cannot use LoopProjectfileProcessor: Loop project file cannot be imported, try installing LoopProjectFile"
    )
from LoopStructural.modelling.features import StructuralFrame
from LoopStructural.modelling.fault import FaultSegment
