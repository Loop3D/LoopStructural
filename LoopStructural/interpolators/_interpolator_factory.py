"""Compatibility shim for the interpolator factory.

The implementation now lives in ``loop_interpolation._interpolator_factory``.
"""

from warnings import warn

from loop_interpolation._interpolator_factory import InterpolatorFactory as _InterpolatorFactory

warn(
    "LoopStructural.interpolators._interpolator_factory is deprecated; use "
    "loop_interpolation._interpolator_factory instead.",
    DeprecationWarning,
    stacklevel=2,
)


class InterpolatorFactory(_InterpolatorFactory):
    """Backward-compatible alias for loop_interpolation.InterpolatorFactory."""
