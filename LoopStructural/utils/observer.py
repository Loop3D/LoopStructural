"""Compatibility re-export: the generic Observer pattern now lives in loop_common."""

from loop_common.observer import Callback, Disposable, Observable, Observer

__all__ = ["Callback", "Disposable", "Observable", "Observer"]
