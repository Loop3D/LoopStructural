"""
Utils
=====
"""

import numpy as np

from ._transformation import EuclideanTransformation
from .exceptions import (
    InterpolatorError,
    LoopException,
    LoopImportError,
    LoopTypeError,
    LoopValueError,
)
from .helper import (
    create_box,
    create_surface,
    get_data_bounding_box,
    get_data_bounding_box_map,
)
from .json_encoder import LoopJSONEncoder
from .logging import get_levels, getLogger, log_to_console, log_to_file

# from ..datatypes._bounding_box import BoundingBox
from .maths import (
    azimuthplunge2vector,
    get_dip_vector,
    get_strike_vector,
    get_vectors,
    normal_vector_to_strike_and_dip,
    rotate,
    strikedip2vector,
)
from .regions import NegativeRegion, PositiveRegion, RegionEverywhere, RegionFunction

rng = np.random.default_rng()

from ._surface import LoopIsosurfacer, surface_list
from .colours import random_colour, random_hex_colour
