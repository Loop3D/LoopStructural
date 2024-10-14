from enum import IntEnum


class SupportType(IntEnum):
    """
    Enum for the different interpolator types

    1-9 should cover interpolators with supports
    9+ are data supported
    """

    StructuredGrid2D = 0
    StructuredGrid = 1
    UnStructuredTetMesh = 2
    P1Unstructured2d = 3
    P2Unstructured2d = 4
    BaseUnstructured2d = 5
    BaseStructured = 6
    TetMesh = 10
    P2UnstructuredTetMesh = 11
    DataSupported = 12


from ._2d_base_unstructured import BaseUnstructured2d
from ._2d_p1_unstructured import P1Unstructured2d
from ._2d_p2_unstructured import P2Unstructured2d
from ._2d_structured_grid import StructuredGrid2D
from ._3d_structured_grid import StructuredGrid
from ._3d_unstructured_tetra import UnStructuredTetMesh
from ._3d_structured_tetra import TetMesh
from ._3d_p2_tetra import P2UnstructuredTetMesh


def no_support(*args, **kwargs):
    return None


support_map = {
    SupportType.StructuredGrid2D: StructuredGrid2D,
    SupportType.StructuredGrid: StructuredGrid,
    SupportType.UnStructuredTetMesh: UnStructuredTetMesh,
    SupportType.P1Unstructured2d: P1Unstructured2d,
    SupportType.P2Unstructured2d: P2Unstructured2d,
    SupportType.TetMesh: TetMesh,
    SupportType.P2UnstructuredTetMesh: P2UnstructuredTetMesh,
    SupportType.DataSupported: no_support,
}

from ._support_factory import SupportFactory

__all__ = [
    "BaseUnstructured2d",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "StructuredGrid2D",
    "StructuredGrid",
    "UnStructuredTetMesh",
    "TetMesh",
    "P2UnstructuredTetMesh",
    "support_map",
    "SupportType",
]
