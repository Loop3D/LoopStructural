"""
Exported file formats
"""

from enum import Enum


class FileFormat(Enum):
    """Enumeration of file export formats"""

    OBJ = 1  # Not supported yet
    VTK = 2
    GOCAD = 3
    GLTF = 4  # Not supported yet
    NUMPY = 5
