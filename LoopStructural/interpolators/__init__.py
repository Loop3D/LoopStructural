"""
Interpolators and interpolation supports

"""
from enum import IntEnum

class InterpolatorType(IntEnum):
    '''
    Enum for the different interpolator types
    
    1-9 should cover interpolators with supports
    9+ are data supported
    '''
    BASE = 0
    BASE_DISCRETE = 1
    FINITE_DIFFERENCE = 2
    DISCRETE_FOLD = 3
    PIECEWISE_LINEAR = 4
    BASE_DATA_SUPPORTED = 10
    SURFE = 11
from LoopStructural.interpolators.geological_interpolator import GeologicalInterpolator
from LoopStructural.interpolators.discrete_interpolator import DiscreteInterpolator
from LoopStructural.interpolators.structured_tetra import TetMesh
from LoopStructural.interpolators.unstructured_tetra import UnStructuredTetMesh
from LoopStructural.interpolators.structured_grid import StructuredGrid
from LoopStructural.interpolators.finite_difference_interpolator import (
    FiniteDifferenceInterpolator,
)
from LoopStructural.interpolators.piecewiselinear_interpolator import (
    PiecewiseLinearInterpolator,
)
from LoopStructural.interpolators.discrete_fold_interpolator import (
    DiscreteFoldInterpolator,
)



