"""
Utils
=====
"""

from loop_common.utils import rng

from ._api_registry import (
    get_registry,
    get_stable_surface,
    public_api,
    register_external_stable,
)
from ._surface import LoopIsosurfacer, surface_list
from ._transformation import EuclideanTransformation
from .colours import random_colour, random_hex_colour
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
from .logging import (
    FileSink,
    LogSink,
    SqliteSink,
    StreamSink,
    add_sink,
    get_levels,
    getLogger,
    log_to_console,
    log_to_file,
    remove_sink,
    timed,
    timed_stage,
)
from .maths import (
    azimuthplunge2vector,
    get_dip_vector,
    get_strike_vector,
    get_vectors,
    normal_vector_to_dip_and_dip_direction,
    normal_vector_to_strike_and_dip,
    plungeazimuth2vector,
    rotate,
    strikedip2vector,
)
from .observer import Callback, Disposable, Observable
from .regions import NegativeRegion, PositiveRegion, RegionEverywhere, RegionFunction
