"""loop-link: send geoscience data to the LoopConnect viewer."""

from __future__ import annotations

__version__ = "0.1.0"

from .errors import (
    BusyError,
    CancelledError,
    InvalidInputError,
    LoopLinkError,
    MissingExtraError,
    MissingFeatureError,
    NoViewerError,
    ObjectTooLargeError,
    Refusal,
    StaleRevisionError,
    TokenError,
    UnsupportedObjectError,
    ViewerClosedError,
    ViewerRefusedError,
)
from .fetch import FetchResult
from .objects import BlockModel, Categorical, Drillhole, Lines, Points, Surface, TetraMesh, Vectors
from .push import PushResult, push_objects
from .selection import Selection, SelectionWatcher
from .viewer import InventoryEntry, Viewer

__all__ = [
    "BlockModel",
    "BusyError",
    "CancelledError",
    "Categorical",
    "Drillhole",
    "FetchResult",
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
    "Selection",
    "SelectionWatcher",
    "StaleRevisionError",
    "Surface",
    "TetraMesh",
    "TokenError",
    "UnsupportedObjectError",
    "Vectors",
    "Viewer",
    "ViewerClosedError",
    "ViewerRefusedError",
    "push_objects",
]
