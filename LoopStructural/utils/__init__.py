"""Compatibility layer for LoopStructural utils."""

from loop_common.utils import (
    getLogger,
    log_to_file,
    log_to_console,
    get_levels,
    LogSink,
    StreamSink,
    FileSink,
    SqliteSink,
    add_sink,
    remove_sink,
    timed_stage,
    timed,
    EuclideanTransformation,
    get_data_bounding_box,
    get_data_bounding_box_map,
    create_surface,
    create_box,
    LoopException,
    LoopImportError,
    InterpolatorError,
    LoopTypeError,
    LoopValueError,
    rng,
)

from .maths import (
    get_dip_vector,
    get_strike_vector,
    get_vectors,
    strikedip2vector,
    plungeazimuth2vector,
    azimuthplunge2vector,
    normal_vector_to_strike_and_dip,
    normal_vector_to_dip_and_dip_direction,
    rotate,
)
from .regions import RegionEverywhere, RegionFunction, NegativeRegion, PositiveRegion
from .json_encoder import LoopJSONEncoder
from ._surface import LoopIsosurfacer, surface_list
from .colours import random_colour, random_hex_colour
from .observer import Callback, Disposable, Observable
from ._api_registry import public_api, get_registry, get_stable_surface
