from loop_common.geometry import (
    BoundingBox,
    PropertyGrid,
    StructuredGrid,
    StructuredGrid2DGeometry,
    StructuredGrid3DGeometry,
    Surface,
    UnstructuredMesh2DGeometry,
    UnstructuredMeshGeometry,
    ValuePoints,
    VectorPoints,
)

from ..utils._api_registry import register_external_stable

for _cls in (BoundingBox, Surface, ValuePoints, VectorPoints):
    register_external_stable(f"LoopStructural.geometry.{_cls.__name__}", _cls.__init__)
del _cls

__all__ = [
    "BoundingBox",
    "PropertyGrid",
    "StructuredGrid",
    "StructuredGrid2DGeometry",
    "StructuredGrid3DGeometry",
    "Surface",
    "UnstructuredMesh2DGeometry",
    "UnstructuredMeshGeometry",
    "ValuePoints",
    "VectorPoints",
]
