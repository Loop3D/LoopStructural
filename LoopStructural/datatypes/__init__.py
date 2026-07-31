"""Deprecated import path.

``BoundingBox``, ``Surface``, ``ValuePoints`` and ``VectorPoints`` moved to
:mod:`LoopStructural.geometry` in v1.6.x. This shim re-exports them so
existing consumers (e.g. the LoopStructural QGIS plugin) keep working, and
will be removed after two minor releases per the versioning policy in
``ROADMAP.md``.
"""

import warnings

from ..geometry import BoundingBox, Surface, ValuePoints, VectorPoints

warnings.warn(
    "LoopStructural.datatypes is deprecated and will be removed in a future "
    "release; import BoundingBox, Surface, ValuePoints and VectorPoints from "
    "LoopStructural.geometry instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "BoundingBox",
    "Surface",
    "ValuePoints",
    "VectorPoints",
]
