"""
Exported file formats
"""
from enum import Enum

class FileFormat(Enum):
    """ Enumeration of file export formats 
    """
    OBJ = 1
    VTK = 2
    GZ = 3 # Not supported yet
    GLTF = 4 # Not supported yet
