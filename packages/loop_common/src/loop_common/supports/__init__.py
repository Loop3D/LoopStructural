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
    P2StructuredTetMesh = 12
    DataSupported = 13
    RectilinearGrid = 14


from ._2d_base_unstructured import BaseUnstructured2d
from ._2d_p1_unstructured import P1Unstructured2d
from ._2d_p2_unstructured import P2Unstructured2d
from ._2d_structured_grid import StructuredGrid2D
from ._3d_p2_tetra import P2UnstructuredTetMesh
from ._3d_rectilinear_grid import RectilinearGrid
from ._3d_structured_grid import StructuredGrid
from ._3d_structured_tetra import TetMesh
from ._3d_unstructured_tetra import UnStructuredTetMesh
from ._p2_structured_tetra import P2TetMesh


def no_support(*args, **kwargs):
    return None


support_map = {
    SupportType.StructuredGrid2D: StructuredGrid2D,
    SupportType.StructuredGrid: StructuredGrid,
    SupportType.RectilinearGrid: RectilinearGrid,
    SupportType.UnStructuredTetMesh: UnStructuredTetMesh,
    SupportType.P1Unstructured2d: P1Unstructured2d,
    SupportType.P2Unstructured2d: P2Unstructured2d,
    SupportType.TetMesh: TetMesh,
    SupportType.P2UnstructuredTetMesh: P2UnstructuredTetMesh,
    SupportType.P2StructuredTetMesh: P2TetMesh,
    SupportType.DataSupported: no_support,
}

from ._support_factory import SupportFactory

__all__ = [
    "BaseUnstructured2d",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "P2UnstructuredTetMesh",
    "RectilinearGrid",
    "StructuredGrid",
    "StructuredGrid2D",
    "SupportType",
    "TetMesh",
    "UnStructuredTetMesh",
    "support_map",
]
