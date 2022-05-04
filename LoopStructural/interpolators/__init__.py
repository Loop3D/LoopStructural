"""
Interpolators and interpolation supports

"""
__all__ = ['InterpolatorType','GeologicalInterpolator','DiscreteInterpolator','FiniteDifferenceInterpolator','PiecewiseLinearInterpolator','DiscreteFoldInterpolator','SurfeRBFInterpolator','P1Interpolator','P2Interpolator','TetMesh',
    'StructuredGrid',
    'UnStructuredTetMesh',
    'P1Unstructured2d',
    'P2Unstructured2d',
    'StructuredGrid2D',
    'P2UnstructuredTetMesh',]
from enum import IntEnum

from LoopStructural.utils import getLogger
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


from LoopStructural.interpolators._geological_interpolator import GeologicalInterpolator
from LoopStructural.interpolators._discrete_interpolator import DiscreteInterpolator
from LoopStructural.interpolators.supports import (
    TetMesh,
    StructuredGrid,
    UnStructuredTetMesh,
    P1Unstructured2d,
    P2Unstructured2d,
    StructuredGrid2D,
    P2UnstructuredTetMesh,
)


from LoopStructural.interpolators._finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from LoopStructural.interpolators.piecewiselinear_interpolator import (
    PiecewiseLinearInterpolator,
)
from LoopStructural.interpolators._discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)

try:
    from LoopStructural.interpolators._surfe_wrapper import SurfeRBFInterpolator
except ImportError:
    logger.warning('Can\'t import surfepy - to install "pip install surfe"')

logger.warning("Using experimental interpolators: P1Interpolator and P2Interpolator")
from ._p1interpolator import P1Interpolator
from ._p2interpolator import P2Interpolator
