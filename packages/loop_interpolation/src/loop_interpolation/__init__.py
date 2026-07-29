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
    "FDFoldInterpolator",
    "SurfeRBFInterpolator",
    "P1Interpolator",
    "P2Interpolator",
    "TetMesh",
    "StructuredGrid",
    "UnStructuredTetMesh",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "StructuredGrid2D",
    "P2UnstructuredTetMesh",
    "ConstantNormP1Interpolator",
    "ConstantNormFDIInterpolator",
    "ConstraintDiagnosticsReport",
    "ConstraintFamilyDiagnostics",
    "RegionCoverageDiagnostics",
    "DirectionalRegularisation",
    "RegularisationConfig",
    "FoldEvent",
    "FourierSeriesFoldRotationAngleProfile",
    "LambdaFoldRotationAngleProfile",
    "FoldRotationType",
    "get_fold_rotation_profile",
]
from ._interpolatortype import InterpolatorType

from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)

from ._geological_interpolator import GeologicalInterpolator
from ._discrete_interpolator import DiscreteInterpolator
from ._diagnostics import (
    ConstraintDiagnosticsReport,
    ConstraintFamilyDiagnostics,
    RegionCoverageDiagnostics,
)
from ._regularisation import DirectionalRegularisation, RegularisationConfig
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


from ._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from ._p1interpolator import (
    P1Interpolator as PiecewiseLinearInterpolator,
)
from ._discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)
from ._fd_fold_interpolator import FDFoldInterpolator
from ._p2interpolator import P2Interpolator
from ._p1interpolator import P1Interpolator
from ._constant_norm import ConstantNormP1Interpolator, ConstantNormFDIInterpolator

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


# Ensure compatibility between the fallback and imported class
SurfeRBFInterpolator = SurfeRBFInterpolator


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

from ._interpolator_factory import InterpolatorFactory
from ._interpolator_builder import InterpolatorBuilder
from ._fold_event import FoldEvent
from .fold_function import (
    FourierSeriesFoldRotationAngleProfile,
    LambdaFoldRotationAngleProfile,
    FoldRotationType,
    get_fold_rotation_profile,
)
