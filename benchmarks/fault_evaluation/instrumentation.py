"""Opt-in call counters and memory measurement for fault evaluation benchmarks.

Nothing in this module changes production code at import time. The counters
patch methods only while an ``EvaluationInstrumentation`` context is active.
"""

# ruff: noqa: ANN401  (the wrapped objects are untyped model and fixture objects)

from __future__ import annotations

import tracemalloc
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


class PeakMemory:
    """Measure the peak memory that Python and NumPy allocate during a call.

    The meter uses ``tracemalloc``. NumPy reports its array memory to it. The
    peak is the highest traced memory above the value at ``start``. Memory that
    the process freed before ``start`` cannot hide an allocation, so the order
    of the cases does not change the result (a process RSS can reuse freed
    pages and show no growth). The meter does not see memory that a C library
    allocates outside Python and NumPy.
    """

    METHOD = "tracemalloc"

    def __init__(self) -> None:
        self._idle = 0
        self._own_trace = False

    def start(self) -> None:
        """Start tracing if needed, and record the present traced memory."""
        if not tracemalloc.is_tracing():
            tracemalloc.start()
            self._own_trace = True
        tracemalloc.reset_peak()
        self._idle = tracemalloc.get_traced_memory()[0]

    def stop(self) -> int:
        """Return the peak extra memory in bytes, and stop a trace that start began."""
        if not tracemalloc.is_tracing():
            return 0
        peak = tracemalloc.get_traced_memory()[1]
        if self._own_trace:
            tracemalloc.stop()
            self._own_trace = False
        return max(peak - self._idle, 0)
