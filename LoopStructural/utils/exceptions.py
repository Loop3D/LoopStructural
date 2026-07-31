"""Compatibility re-export: LoopStructural's exception hierarchy now lives in loop_common."""

from loop_common.utils import (
    InterpolatorError,
    LoopException,
    LoopImportError,
    LoopTypeError,
    LoopValueError,
)

__all__ = [
    "InterpolatorError",
    "LoopException",
    "LoopImportError",
    "LoopTypeError",
    "LoopValueError",
]
