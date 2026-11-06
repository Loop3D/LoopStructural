"""loop-link: send geoscience data to the LoopConnect viewer."""

from __future__ import annotations

__version__ = "0.1.0"

from .errors import (
    BusyError,
    InvalidInputError,
    LoopLinkError,
    MissingExtraError,
    MissingFeatureError,
    NoViewerError,
    ObjectTooLargeError,
    Refusal,
    TokenError,
    ViewerClosedError,
    ViewerRefusedError,
)
from .objects import BlockModel, Categorical, Drillhole, Lines, Points, Surface, TetraMesh, Vectors
from .push import PushResult, push_objects
from .viewer import InventoryEntry, Viewer

__all__ = [
    "BlockModel",
    "BusyError",
    "Categorical",
    "Drillhole",
    "InvalidInputError",
    "InventoryEntry",
    "Lines",
    "LoopLinkError",
    "MissingExtraError",
    "MissingFeatureError",
    "NoViewerError",
    "ObjectTooLargeError",
    "Points",
    "PushResult",
    "Refusal",
    "Surface",
    "TetraMesh",
    "TokenError",
    "Vectors",
    "Viewer",
    "ViewerClosedError",
    "ViewerRefusedError",
    "push_objects",
]
