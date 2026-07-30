from loop_common.geometry import (
    Surface,
    ValuePoints,
    VectorPoints,
    StructuredGrid3DGeometry,
    StructuredGrid2DGeometry,
    UnstructuredMeshGeometry,
    UnstructuredMesh2DGeometry,
)

# BoundingBox is kept as the local implementation (global_origin/global_maximum
# reprojection API) rather than loop_common's (local_origin/local_rotation
# affine-transform API) -- the two have diverged onto incompatible constructor
# signatures and property names (see ROADMAP.md 2c-3) and core callers like
# GeologicalModel still depend on the global_origin/global_maximum surface.
from ._bounding_box import BoundingBox
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
