"""Compatibility re-export: the generic Observer pattern now lives in loop_common."""

from loop_common.observer import Callback, Disposable, Observable, Observer

from ._api_registry import register_external_stable

__all__ = ["Callback", "Disposable", "Observable", "Observer"]

register_external_stable("LoopStructural.utils.observer.Observable", Observable.__init__)
