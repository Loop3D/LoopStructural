"""
Geological modelling classes and functions

"""

__all__ = [
    "GeologicalModel",
    "ProcessInputData",
    "Map2LoopProcessor",
    "LoopProjectfileProcessor",
]
from ..utils import LoopImportError, getLogger
from .core.geological_model import GeologicalModel

logger = getLogger(__name__)
from ..modelling.input import (
    Map2LoopProcessor,
    ProcessInputData,
)

try:
    from ..modelling.input.project_file import LoopProjectfileProcessor
except (LoopImportError, ImportError):
    logger.warning(
        "Cannot use LoopProjectfileProcessor: Loop project file cannot be imported, try installing LoopProjectFile"
    )
