"""Interpolators and interpolation supports for LoopStructural.

This module provides various interpolation methods and support structures
for geological modelling, including finite difference, piecewise linear,
and radial basis function interpolators.
"""

__all__ = [
    "DiscreteFoldInterpolator",
    "DiscreteInterpolator",
    "FiniteDifferenceInterpolator",
    "GeologicalInterpolator",
    "InterpolatorBuilder",
    "InterpolatorFactory",
    "InterpolatorType",
    "Operator",
    "P1Interpolator",
    "P1Unstructured2d",
    "P2Interpolator",
    "P2Unstructured2d",
    "P2UnstructuredTetMesh",
    "PiecewiseLinearInterpolator",
    "StructuredGrid",
    "StructuredGrid2D",
    "StructuredGridSupport",
    "SupportType",
    "SurfeRBFInterpolator",
    "TetMesh",
    "UnStructuredTetMesh",
    "add_element_anisotropy_constraints",
    "add_fd_anisotropy_constraints",
]
from loop_common.supports import (
    P1Unstructured2d,
    P2Unstructured2d,
    P2UnstructuredTetMesh,
    StructuredGrid,
    StructuredGrid2D,
    SupportType,
    TetMesh,
    UnStructuredTetMesh,
)
from loop_interpolation import (
    ConstantNormFDIInterpolator,
    ConstantNormP1Interpolator,
    DiscreteInterpolator,
    FiniteDifferenceInterpolator,
    GeologicalInterpolator,
    InterpolatorBuilder,
    InterpolatorFactory,
    InterpolatorType,
    P1Interpolator,
    P2Interpolator,
    PiecewiseLinearInterpolator,
    SurfeRBFInterpolator,
    add_element_anisotropy_constraints,
    add_fd_anisotropy_constraints,
    interpolator_map,
    interpolator_string_map,
    support_interpolator_map,
)
from loop_interpolation._operator import Operator

from ..utils import getLogger

logger = getLogger(__name__)

# Legacy LoopStructural name kept for backwards compatibility.
StructuredGridSupport = StructuredGrid


def __getattr__(name):
    # Deferred to avoid a LoopStructural.interpolators <-> LoopStructural.modelling
    # import cycle: DiscreteFoldInterpolator (the fold-aware P1Interpolator
    # subclass) lives in modelling.features.fold, which itself depends on
    # modelling being importable. Resolving it lazily here means this module
    # can still be imported first, as LoopStructural/__init__.py does.
    if name == "DiscreteFoldInterpolator":
        from ..modelling.features.fold import DiscreteFoldInterpolator

        return DiscreteFoldInterpolator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
