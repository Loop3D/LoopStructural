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
    DiscreteFoldInterpolator,
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
    interpolator_map,
    interpolator_string_map,
    support_interpolator_map,
)
from loop_interpolation._operator import Operator

from ..utils import getLogger

logger = getLogger(__name__)

# Legacy LoopStructural name kept for backwards compatibility.
StructuredGridSupport = StructuredGrid
