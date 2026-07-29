"""Compatibility shim for finite-difference interpolation.

The implementation now lives in ``loop_interpolation._finite_difference_interpolator``.
"""

from warnings import warn

from loop_interpolation._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
    compute_weighting,
)

warn(
    "LoopStructural.interpolators._finite_difference_interpolator is deprecated; use "
    "loop_interpolation._finite_difference_interpolator instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["FiniteDifferenceInterpolator", "compute_weighting"]
