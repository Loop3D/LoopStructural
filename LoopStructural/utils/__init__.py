"""
Utils
=====
"""

from .logging import getLogger, log_to_file, log_to_console
from .exceptions import (
    LoopException,
    LoopImportError,
    InterpolatorError,
    LoopTypeError,
    LoopValueError,
)
from ._transformation import EuclideanTransformation
from .helper import (
    get_data_bounding_box,
    get_data_bounding_box_map,
)

# from ..datatypes._bounding_box import BoundingBox
from .maths import (
    get_dip_vector,
    get_strike_vector,
    get_vectors,
    strikedip2vector,
    azimuthplunge2vector,
    normal_vector_to_strike_and_dip,
    rotate,
)
from .helper import create_surface, create_box
from .regions import RegionEverywhere, RegionFunction, NegativeRegion, PositiveRegion

from .json_encoder import LoopJSONEncoder
import numpy as np

rng = np.random.default_rng()
