"""
Interpolators and interpolation supports

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
    "TetMesh",
    "StructuredGrid",
    "UnStructuredTetMesh",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "StructuredGrid2D",
    "P2UnstructuredTetMesh",
]
from enum import IntEnum

from ..utils import getLogger

logger = getLogger(__name__)


class InterpolatorType(IntEnum):
    """
    Enum for the different interpolator types

    1-9 should cover interpolators with supports
    9+ are data supported
    """

    BASE = 0
    BASE_DISCRETE = 1
    FINITE_DIFFERENCE = 2
    DISCRETE_FOLD = 3
    PIECEWISE_LINEAR = 4
    PIECEWISE_QUADRATIC = 5
    BASE_DATA_SUPPORTED = 10
    SURFE = 11


interpolator_string_map = {
    "FDI": InterpolatorType.FINITE_DIFFERENCE,
    "PLI": InterpolatorType.PIECEWISE_LINEAR,
    "P2": InterpolatorType.PIECEWISE_QUADRATIC,
    "P1": InterpolatorType.PIECEWISE_LINEAR,
    "DFI": InterpolatorType.DISCRETE_FOLD,
    'surfe': InterpolatorType.SURFE,
}
from ..interpolators._geological_interpolator import GeologicalInterpolator
from ..interpolators._discrete_interpolator import DiscreteInterpolator
from ..interpolators.supports import (
    TetMesh,
    StructuredGrid,
    UnStructuredTetMesh,
    P1Unstructured2d,
    P2Unstructured2d,
    StructuredGrid2D,
    P2UnstructuredTetMesh,
    SupportType,
)


from ..interpolators._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from ..interpolators._p1interpolator import (
    P1Interpolator as PiecewiseLinearInterpolator,
)
from ..interpolators._discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)
from ..interpolators._p2interpolator import P2Interpolator
from ..interpolators._p1interpolator import P1Interpolator

try:
    from ..interpolators._surfe_wrapper import SurfeRBFInterpolator
except ImportError:
    logger.warning("Surfe is not installed, SurfeRBFInterpolator will not be available")
    SurfeRBFInterpolator = None
interpolator_map = {
    InterpolatorType.BASE: GeologicalInterpolator,
    InterpolatorType.BASE_DISCRETE: DiscreteInterpolator,
    InterpolatorType.FINITE_DIFFERENCE: FiniteDifferenceInterpolator,
    InterpolatorType.DISCRETE_FOLD: DiscreteFoldInterpolator,
    InterpolatorType.PIECEWISE_LINEAR: P1Interpolator,
    InterpolatorType.PIECEWISE_QUADRATIC: P2Interpolator,
    InterpolatorType.BASE_DATA_SUPPORTED: GeologicalInterpolator,
    InterpolatorType.SURFE: SurfeRBFInterpolator,
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
}

from ._interpolator_factory import InterpolatorFactory
from ._interpolator_builder import InterpolatorBuilder
# from ._api import LoopInterpolator
