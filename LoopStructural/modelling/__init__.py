"""
Geological modelling classes and functions

"""

__all__ = [
    "GeologicalModel",
    "ProcessInputData",
    "Map2LoopProcessor",
    "LoopProjectfileProcessor",
]
from ..utils import getLogger
from ..utils import LoopImportError
from .core.geological_model import GeologicalModel

logger = getLogger(__name__)
from ..modelling.input import (
    ProcessInputData,
    Map2LoopProcessor,
)

try:
    from ..modelling.input.project_file import LoopProjectfileProcessor
except (LoopImportError, ImportError):
    class LoopProjectfileProcessor(ProcessInputData):
        """
        Dummy class to handle the case where LoopProjectFile is not installed.
        This will raise a warning when used.
        """

        def __init__(self, *args, **kwargs):
            raise LoopImportError(
                "LoopProjectFile cannot be imported. Please install LoopProjectFile."
            )
    
