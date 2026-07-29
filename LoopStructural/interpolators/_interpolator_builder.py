"""Compatibility shim for the fluent interpolator builder.

The implementation now lives in ``loop_interpolation._interpolator_builder``.
"""

from warnings import warn

from loop_interpolation._interpolator_builder import InterpolatorBuilder as _InterpolatorBuilder

warn(
    "LoopStructural.interpolators._interpolator_builder is deprecated; use "
    "loop_interpolation._interpolator_builder instead. This compatibility shim "
    "will be removed in 2 minor releases.",
    DeprecationWarning,
    stacklevel=2,
)


class InterpolatorBuilder(_InterpolatorBuilder):
    """Backward-compatible alias for loop_interpolation.InterpolatorBuilder."""
