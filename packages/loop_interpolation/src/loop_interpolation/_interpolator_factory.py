"""Interpolator creation and reconstruction backend.

InterpolatorFactory is the single place that knows how to:
- map type identifiers to concrete interpolator classes,
- create supports from bounding boxes, and
- reconstruct interpolators from serialized dictionaries.

Higher-level fluent APIs should delegate construction to this module.
"""

from typing import Optional, Union

from loop_common.geometry import BoundingBox
from loop_common.supports import SupportFactory

from . import (
    interpolator_map,
    InterpolatorType,
    support_interpolator_map,
    interpolator_string_map,
)
import numpy as np


class InterpolatorFactory:
    """Authoritative constructor/reconstructor for interpolation objects."""

    @staticmethod
    def _normalise_interpolator_type(
        interpolatortype: Union[str, InterpolatorType],
    ) -> InterpolatorType:
        if isinstance(interpolatortype, str):
            if interpolatortype in interpolator_string_map:
                return interpolator_string_map[interpolatortype]
            if interpolatortype in InterpolatorType.__members__:
                return InterpolatorType[interpolatortype]
            return InterpolatorType(interpolatortype)
        return interpolatortype

    @staticmethod
    def create_interpolator(
        interpolatortype: Optional[Union[str, InterpolatorType]] = None,
        boundingbox: Optional[BoundingBox] = None,
        nelements: Optional[int] = None,
        element_volume: Optional[float] = None,
        support=None,
        buffer: Optional[float] = None,
        solver: Optional[str] = None,
    ):
        if interpolatortype is None:
            raise ValueError("No interpolator type specified")
        if support is None and boundingbox is None:
            raise ValueError("No bounding box specified")

        interpolatortype = InterpolatorFactory._normalise_interpolator_type(interpolatortype)
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
        interpolator = interpolator_map[interpolatortype](support)
        if boundingbox is not None:
            interpolator.bounding_box = boundingbox
        if solver is not None:
            interpolator.solver = solver
        return interpolator

    @staticmethod
    def from_dict(d):
        d = d.copy()
        interpolator_type = d.pop("type", None)
        if interpolator_type is None:
            raise ValueError("No interpolator type specified")
        interpolator_type = InterpolatorFactory._normalise_interpolator_type(interpolator_type)

        support = d.pop("support", None)
        if isinstance(support, dict):
            support_payload = support.copy()
            try:
                support = SupportFactory.from_dict(support_payload)
            except TypeError as exc:
                # Some support classes (e.g., TetMesh) do not accept rotation_xy in __init__.
                if "rotation_xy" in support_payload and "rotation_xy" in str(exc):
                    support_payload.pop("rotation_xy", None)
                    support = SupportFactory.from_dict(support_payload)
                else:
                    raise

        data = d.pop("data", None)
        c = d.pop("c", None)
        up_to_date = bool(d.pop("up_to_date", False))
        valid = bool(d.pop("valid", True))

        interpolator = InterpolatorFactory.create_interpolator(
            interpolator_type,
            support=support,
            **d,
        )

        if data is not None:
            if data.get("value") is not None:
                arr = np.asarray(data["value"], dtype=float)
                if arr.size > 0:
                    interpolator.set_value_constraints(arr)
            if data.get("gradient") is not None:
                arr = np.asarray(data["gradient"], dtype=float)
                if arr.size > 0:
                    interpolator.set_gradient_constraints(arr)
            if data.get("normal") is not None:
                arr = np.asarray(data["normal"], dtype=float)
                if arr.size > 0:
                    interpolator.set_normal_constraints(arr)
            if data.get("tangent") is not None:
                arr = np.asarray(data["tangent"], dtype=float)
                if arr.size > 0:
                    interpolator.set_tangent_constraints(arr)
            if data.get("interface") is not None:
                arr = np.asarray(data["interface"], dtype=float)
                if arr.size > 0:
                    interpolator.set_interface_constraints(arr)
            if data.get("inequality") is not None:
                arr = np.asarray(data["inequality"], dtype=float)
                if arr.size > 0:
                    interpolator.set_value_inequality_constraints(arr)
            if data.get("inequality_pairs") is not None:
                arr = np.asarray(data["inequality_pairs"], dtype=float)
                if arr.size > 0:
                    interpolator.set_inequality_pairs_constraints(arr)

        if c is not None and hasattr(interpolator, "c"):
            interpolator.c = np.asarray(c, dtype=float)

        interpolator.up_to_date = up_to_date
        interpolator.valid = valid
        return interpolator

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
