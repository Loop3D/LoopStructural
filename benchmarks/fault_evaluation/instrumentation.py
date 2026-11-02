"""Opt-in call counters and memory sampling for fault evaluation benchmarks.

Nothing in this module changes production code at import time. The counters
patch methods only while an ``EvaluationInstrumentation`` context is active.
"""

# ruff: noqa: ANN401  (the wrapped objects are untyped model and fixture objects)

from __future__ import annotations

import sys
import threading
from typing import Any, Callable

from loop_interpolation._discrete_interpolator import DiscreteInterpolator

from LoopStructural.modelling.features.fault import FaultSegment

_RSS_BYTES_PER_UNIT = {"bytes": 1, "kibibytes": 1024}


def rss_to_bytes(value: float, unit: str) -> int:
    """Convert a native memory value to bytes. Reject an unknown unit."""
    try:
        factor = _RSS_BYTES_PER_UNIT[unit]
    except KeyError:
        raise ValueError(
            f"Unknown RSS unit {unit!r}; expected one of {sorted(_RSS_BYTES_PER_UNIT)}"
        ) from None
    return int(value * factor)


class EvaluationInstrumentation:
    """Count fault restorations and local field evaluations inside a context."""

    def __init__(self) -> None:
        self._saved: list[tuple[type, str, Any]] = []
        self._depth = 0
        self.reset()

    def reset(self) -> None:
        """Set all counters to zero."""
        self._counts = {
            "fault_calls_total": 0,
            "fault_calls_top_level": 0,
            "fault_points_total": 0,
            "field_calls_total": 0,
            "field_points_total": 0,
        }

    def snapshot(self) -> dict[str, int]:
        """Return a copy of the counters."""
        return dict(self._counts)

    def _wrap_fault(self, original: Callable) -> Callable:
        def apply_to_points(fault: Any, points: Any, *args: Any, **kwargs: Any) -> Any:
            self._counts["fault_calls_total"] += 1
            self._counts["fault_points_total"] += len(points)
            if self._depth == 0:
                self._counts["fault_calls_top_level"] += 1
            self._depth += 1
            try:
                return original(fault, points, *args, **kwargs)
            finally:
                self._depth -= 1

        return apply_to_points

    def _wrap_field(self, original: Callable) -> Callable:
        def evaluate_local(interpolator: Any, locations: Any, *args: Any, **kwargs: Any) -> Any:
            self._counts["field_calls_total"] += 1
            self._counts["field_points_total"] += len(locations)
            return original(interpolator, locations, *args, **kwargs)

        return evaluate_local

    def _patch(self, owner: type, name: str, wrapper_factory: Callable) -> None:
        original = owner.__dict__[name]
        self._saved.append((owner, name, original))
        setattr(owner, name, wrapper_factory(original))

    def __enter__(self) -> EvaluationInstrumentation:  # noqa: PYI034
        """Patch the evaluation methods."""
        if self._saved:
            raise RuntimeError("Instrumentation is already active")
        try:
            self._patch(FaultSegment, "apply_to_points", self._wrap_fault)
            self._patch(DiscreteInterpolator, "_evaluate_value_local", self._wrap_field)
            self._patch(DiscreteInterpolator, "_evaluate_gradient_local", self._wrap_field)
        except BaseException:
            self._restore()
            raise
        return self

    def _restore(self) -> None:
        while self._saved:
            owner, name, original = self._saved.pop()
            setattr(owner, name, original)
        self._depth = 0

    def __exit__(self, *exc_info: object) -> None:
        """Restore the original methods, also after an error."""
        self._restore()


def _current_rss_bytes() -> int:
    """Return the resident memory of this process in bytes."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class _MemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        kernel32 = ctypes.WinDLL("kernel32")
        psapi = ctypes.WinDLL("psapi")
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [
            wintypes.HANDLE,
            ctypes.POINTER(_MemoryCounters),
            wintypes.DWORD,
        ]
        counters = _MemoryCounters()
        counters.cb = ctypes.sizeof(_MemoryCounters)
        if not psapi.GetProcessMemoryInfo(
            kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
        ):
            raise OSError("GetProcessMemoryInfo failed")
        return rss_to_bytes(counters.WorkingSetSize, "bytes")
    try:
        import os

        with open("/proc/self/statm", encoding="utf8") as handle:
            pages = int(handle.read().split()[1])
        return rss_to_bytes(pages * os.sysconf("SC_PAGE_SIZE"), "bytes")
    except (OSError, ValueError, IndexError):
        import resource

        # Without /proc only the process high-water mark is available.
        unit = "bytes" if sys.platform == "darwin" else "kibibytes"
        return rss_to_bytes(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, unit)


class RssSampler:
    """Sample process memory in a thread and report the peak above the idle value."""

    def __init__(self, interval_seconds: float = 0.002) -> None:
        self._interval = interval_seconds
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._idle = 0
        self._peak = 0

    def _sample(self) -> None:
        self._peak = max(self._peak, _current_rss_bytes())

    def _run(self) -> None:
        while not self._stop.wait(self._interval):
            self._sample()

    def start(self) -> None:
        """Record the idle value and start sampling."""
        self._stop.clear()
        self._idle = self._peak = _current_rss_bytes()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> int:
        """Stop sampling. Return the peak extra memory in bytes."""
        if self._thread is None:
            return 0
        self._stop.set()
        self._thread.join()
        self._thread = None
        self._sample()
        return max(self._peak - self._idle, 0)
