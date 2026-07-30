from loop_common.geometry import (
    BoundingBox,
    Surface,
    ValuePoints,
    VectorPoints,
    StructuredGrid3DGeometry,
    StructuredGrid2DGeometry,
    UnstructuredMeshGeometry,
    UnstructuredMesh2DGeometry,
)

from ._structured_grid import StructuredGrid

from ..utils._api_registry import register_external_stable

for _cls in (BoundingBox, Surface, ValuePoints, VectorPoints):
    register_external_stable(f"LoopStructural.geometry.{_cls.__name__}", _cls.__init__)
del _cls

__all__ = [
    "Surface",
    "BoundingBox",
    "ValuePoints",
    "VectorPoints",
    "StructuredGrid",
    "StructuredGrid3DGeometry",
    "StructuredGrid2DGeometry",
    "UnstructuredMeshGeometry",
    "UnstructuredMesh2DGeometry",
]
