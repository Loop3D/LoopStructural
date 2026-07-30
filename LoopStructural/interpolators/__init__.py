"""Interpolators and interpolation supports for LoopStructural.

This module provides various interpolation methods and support structures
for geological modelling, including finite difference, piecewise linear,
and radial basis function interpolators.
"""

__all__ = [
    "InterpolatorType",
    "GeologicalInterpolator",
    "DiscreteInterpolator",
    "FiniteDifferenceInterpolator",
    "PiecewiseLinearInterpolator",
    "DiscreteFoldInterpolator",
    "SurfeRBFInterpolator",
    "P1Interpolator",
    "P2Interpolator",
    "Operator",
    "TetMesh",
    "StructuredGridSupport",
    "StructuredGrid",
    "UnStructuredTetMesh",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "StructuredGrid2D",
    "P2UnstructuredTetMesh",
    "SupportType",
    "InterpolatorFactory",
    "InterpolatorBuilder",
]
from loop_interpolation import (
    InterpolatorType,
    GeologicalInterpolator,
    DiscreteInterpolator,
    FiniteDifferenceInterpolator,
    PiecewiseLinearInterpolator,
    DiscreteFoldInterpolator,
    SurfeRBFInterpolator,
    P1Interpolator,
    P2Interpolator,
    ConstantNormP1Interpolator,
    ConstantNormFDIInterpolator,
    InterpolatorFactory,
    InterpolatorBuilder,
    interpolator_map,
    interpolator_string_map,
    support_interpolator_map,
)
from loop_common.supports import (
    TetMesh,
    StructuredGrid,
    UnStructuredTetMesh,
    P1Unstructured2d,
    P2Unstructured2d,
    StructuredGrid2D,
    P2UnstructuredTetMesh,
    SupportType,
)
from loop_interpolation._operator import Operator

from ..utils import getLogger

logger = getLogger(__name__)

# Legacy LoopStructural name kept for backwards compatibility.
StructuredGridSupport = StructuredGrid
