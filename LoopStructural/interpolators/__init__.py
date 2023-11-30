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
import LoopStructural

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
)


from ..interpolators._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from ..interpolators.piecewiselinear_interpolator import (
    PiecewiseLinearInterpolator,
)
from ..interpolators._discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)

try:
    from ..interpolators._surfe_wrapper import SurfeRBFInterpolator
except ImportError:
    logger.warning('Can\'t import surfepy - to install "pip install surfe"')

logger.warning("Using experimental interpolators: P1Interpolator and P2Interpolator")
from ._p1interpolator import P1Interpolator
from ._p2interpolator import P2Interpolator
from ._builders import get_interpolator

interpolator_map = {
    InterpolatorType.BASE: GeologicalInterpolator,
    InterpolatorType.BASE_DISCRETE: DiscreteInterpolator,
    InterpolatorType.FINITE_DIFFERENCE: FiniteDifferenceInterpolator,
    InterpolatorType.DISCRETE_FOLD: DiscreteFoldInterpolator,
    InterpolatorType.PIECEWISE_LINEAR: PiecewiseLinearInterpolator,
    InterpolatorType.PIECEWISE_QUADRATIC: PiecewiseLinearInterpolator,
    InterpolatorType.BASE_DATA_SUPPORTED: GeologicalInterpolator,
    # InterpolatorType.SURFE: SurfeRBFInterpolator,
}
from ._interpolator_factory import InterpolatorFactory
