from typing import Optional, Union
from .supports import SupportFactory
from . import (
    interpolator_map,
    InterpolatorType,
    support_interpolator_map,
    interpolator_string_map,
)
from LoopStructural.datatypes import BoundingBox
import numpy as np


class InterpolatorFactory:
    @staticmethod
    def create_interpolator(
        interpolatortype: Optional[Union[str, InterpolatorType]] = None,
        boundingbox: Optional[BoundingBox] = None,
        nelements: Optional[int] = None,
        element_volume: Optional[float] = None,
        support=None,
        buffer: Optional[float] = None,
    ):
        if interpolatortype is None:
            raise ValueError("No interpolator type specified")
        if boundingbox is None:
            raise ValueError("No bounding box specified")

        if isinstance(interpolatortype, str):
            interpolatortype = interpolator_string_map[interpolatortype]
        if support is None:
            # raise Exception("Support must be specified")

            supporttype = support_interpolator_map[interpolatortype][boundingbox.dimensions]

            support = SupportFactory.create_support_from_bbox(
                supporttype,
                bounding_box=boundingbox,
                nelements=nelements,
                element_volume=element_volume,
                buffer=buffer,
            )
        return interpolator_map[interpolatortype](support)

    @staticmethod
    def from_dict(d):
        d = d.copy()
        interpolator_type = d.pop("type", None)
        if interpolator_type is None:
            raise ValueError("No interpolator type specified")
        return InterpolatorFactory.create_interpolator(interpolator_type, **d)

    @staticmethod
    def get_supported_interpolators():
        return interpolator_map.keys()

    @staticmethod
    def create_interpolator_with_data(
        interpolatortype: str,
        boundingbox: BoundingBox,
        nelements: int,
        element_volume: Optional[float] = None,
        support=None,
        value_constraints: Optional[np.ndarray] = None,
        gradient_norm_constraints: Optional[np.ndarray] = None,
        gradient_constraints: Optional[np.ndarray] = None,
    ):
        interpolator = InterpolatorFactory.create_interpolator(
            interpolatortype, boundingbox, nelements, element_volume, support
        )
        if value_constraints is not None:
            interpolator.set_value_constraints(value_constraints)
        if gradient_norm_constraints is not None:
            interpolator.set_normal_constraints(gradient_norm_constraints)
        if gradient_constraints is not None:
            interpolator.set_gradient_constraints(gradient_constraints)
        interpolator.setup()
        return interpolator
