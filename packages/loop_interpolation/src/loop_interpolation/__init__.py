"""Interpolators and interpolation supports for LoopStructural.

This module provides various interpolation methods and support structures
for geological modelling, including finite difference, piecewise linear,
and radial basis function interpolators.
"""

__all__ = [
    "ConstantNormFDIInterpolator",
    "ConstantNormP1Interpolator",
    "ConstraintDiagnosticsReport",
    "ConstraintFamilyDiagnostics",
    "DirectionalRegularisation",
    "DiscreteFoldInterpolator",
    "DiscreteInterpolator",
    "FDFoldInterpolator",
    "FiniteDifferenceInterpolator",
    "FoldEvent",
    "FoldRotationType",
    "FourierSeriesFoldRotationAngleProfile",
    "GeologicalInterpolator",
    "InterpolatorType",
    "LambdaFoldRotationAngleProfile",
    "P1Interpolator",
    "P1Unstructured2d",
    "P2Interpolator",
    "P2Unstructured2d",
    "P2UnstructuredTetMesh",
    "PiecewiseLinearInterpolator",
    "RegionCoverageDiagnostics",
    "RegularisationConfig",
    "StructuredGrid",
    "StructuredGrid2D",
    "SurfeRBFInterpolator",
    "TetMesh",
    "UnStructuredTetMesh",
    "get_fold_rotation_profile",
]
from loop_common.logging import get_logger as getLogger

from ._interpolatortype import InterpolatorType

logger = getLogger(__name__)

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

from ._constant_norm import ConstantNormFDIInterpolator, ConstantNormP1Interpolator
from ._diagnostics import (
    ConstraintDiagnosticsReport,
    ConstraintFamilyDiagnostics,
    RegionCoverageDiagnostics,
)
from ._discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)
from ._discrete_interpolator import DiscreteInterpolator
from ._fd_fold_interpolator import FDFoldInterpolator
from ._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from ._geological_interpolator import GeologicalInterpolator
from ._p1interpolator import P1Interpolator
from ._p1interpolator import (
    P1Interpolator as PiecewiseLinearInterpolator,
)
from ._p2interpolator import P2Interpolator
from ._regularisation import DirectionalRegularisation, RegularisationConfig

try:
    from ._surfe_wrapper import SurfeRBFInterpolator
except ImportError:

    class SurfeRBFInterpolator(GeologicalInterpolator):
        """
        Dummy class to handle the case where Surfe is not installed.
        This will raise a warning when used.
        """

        def __new__(cls, *args, **kwargs):
            raise ImportError(
                "Surfe cannot be imported. Please install Surfe. pip install surfe/ conda install -c loop3d surfe"
            )


interpolator_string_map = {
    "FDI": InterpolatorType.FINITE_DIFFERENCE,
    "PLI": InterpolatorType.PIECEWISE_LINEAR,
    "P2": InterpolatorType.PIECEWISE_QUADRATIC,
    "P1": InterpolatorType.PIECEWISE_LINEAR,
    "DFI": InterpolatorType.DISCRETE_FOLD,
    "surfe": InterpolatorType.SURFE,
    "FDI_CN": InterpolatorType.FINITE_DIFFERENCE_CONSTANT_NORM,
    "P1_CN": InterpolatorType.PIECEWISE_LINEAR_CONSTANT_NORM,
}

# Define the mapping after all imports
interpolator_map = {
    InterpolatorType.BASE: GeologicalInterpolator,
    InterpolatorType.BASE_DISCRETE: DiscreteInterpolator,
    InterpolatorType.FINITE_DIFFERENCE: FiniteDifferenceInterpolator,
    InterpolatorType.DISCRETE_FOLD: DiscreteFoldInterpolator,
    InterpolatorType.PIECEWISE_LINEAR: P1Interpolator,
    InterpolatorType.PIECEWISE_QUADRATIC: P2Interpolator,
    InterpolatorType.BASE_DATA_SUPPORTED: GeologicalInterpolator,
    InterpolatorType.SURFE: SurfeRBFInterpolator,
    InterpolatorType.PIECEWISE_LINEAR_CONSTANT_NORM: ConstantNormP1Interpolator,
    InterpolatorType.FINITE_DIFFERENCE_CONSTANT_NORM: ConstantNormFDIInterpolator,
}

support_interpolator_map = {
    InterpolatorType.FINITE_DIFFERENCE: {
        2: SupportType.StructuredGrid2D,
        3: SupportType.StructuredGrid,
    },
    InterpolatorType.DISCRETE_FOLD: {3: SupportType.TetMesh, 2: SupportType.P1Unstructured2d},
    InterpolatorType.PIECEWISE_LINEAR: {3: SupportType.TetMesh, 2: SupportType.P1Unstructured2d},
    InterpolatorType.PIECEWISE_QUADRATIC: {
        3: SupportType.P2UnstructuredTetMesh,
        2: SupportType.P2Unstructured2d,
    },
    InterpolatorType.SURFE: {
        3: SupportType.DataSupported,
        2: SupportType.DataSupported,
    },
    InterpolatorType.PIECEWISE_LINEAR_CONSTANT_NORM: {
        3: SupportType.TetMesh,
        2: SupportType.P1Unstructured2d,
    },
    InterpolatorType.FINITE_DIFFERENCE_CONSTANT_NORM: {
        3: SupportType.StructuredGrid,
        2: SupportType.StructuredGrid2D,
    },
}

from ._fold_event import FoldEvent
from ._interpolator_builder import InterpolatorBuilder
from ._interpolator_factory import InterpolatorFactory
from .fold_function import (
    FoldRotationType,
    FourierSeriesFoldRotationAngleProfile,
    LambdaFoldRotationAngleProfile,
    get_fold_rotation_profile,
)
