"""Interpolator creation and reconstruction backend.

InterpolatorFactory is the single place that knows how to:
- map type identifiers to concrete interpolator classes,
- create supports from bounding boxes, and
- reconstruct interpolators from serialized dictionaries.

Higher-level fluent APIs should delegate construction to this module.
"""

from __future__ import annotations

import numpy as np
from loop_common.geometry import BoundingBox
from loop_common.supports import SupportFactory

from . import (
    InterpolatorType,
    interpolator_map,
    interpolator_string_map,
    support_interpolator_map,
)
from ._frame_persistence import InterpolatorFormatError, payload_kind, resolve_frame
from ._validation import ValidationError


class InterpolatorFactory:
    """Authoritative constructor/reconstructor for interpolation objects."""

    @staticmethod
    def _normalise_interpolator_type(
        interpolatortype: str | InterpolatorType,
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
        interpolatortype: str | InterpolatorType | None = None,
        boundingbox: BoundingBox | None = None,
        nelements: int | None = None,
        element_volume: float | None = None,
        support=None,
        buffer: float | None = None,
        solver: str | None = None,
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
    def from_dict(d, *, local_origin=None, local_rotation=None, declare_identity=False):
        """Build an interpolator from a dict.

        The function accepts two kinds of dict. A dict is a saved payload if
        and only if ``d.get("support")`` is a dict. Any other dict is a
        creation call and keeps the old code path with no change.

        For a saved payload the function checks the kind of the payload first,
        then builds the support and the interpolator, resolves the frame of the
        payload, and restores the stored constraints without a conversion.

        Parameters
        ----------
        d : dict
            Creation call or saved payload. The function does not change it.
        local_origin : array-like, optional
            Origin of the local frame, given by the caller of a saved payload.
        local_rotation : array-like, optional
            Rotation of the local frame, given by the caller of a saved payload.
        declare_identity : bool, default False
            If True, the caller declares an identity frame for a saved payload.

        Returns
        -------
        GeologicalInterpolator
            The new interpolator. If a step raises, no interpolator is returned.

        Raises
        ------
        ValueError
            If the creation call has no interpolator type.
        InterpolatorFormatError
            If a frame argument is given with a creation call, if the saved
            payload is malformed, or if its frame cannot be resolved.
        """
        caller_frame = (
            local_origin is not None or local_rotation is not None or bool(declare_identity)
        )
        saved = isinstance(d.get("support"), dict)
        if not saved and caller_frame:
            raise InterpolatorFormatError(
                "malformed", "frame arguments apply only to a saved payload"
            )
        kind = payload_kind(d, caller_frame) if saved else None
        d = d.copy()
        frame_keys = {}
        if saved:
            for key in ("format_version", "constraint_space", "transform", "bounding_box"):
                if key in d:
                    frame_keys[key] = d.pop(key)
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

        if saved:
            interpolator.bounding_box = resolve_frame(
                frame_keys,
                kind,
                support,
                interpolator.dimensions,
                local_origin=local_origin,
                local_rotation=local_rotation,
                declare_identity=declare_identity,
            )
            if data is not None:
                try:
                    interpolator._restore_local_constraints(data)
                except (ValueError, ValidationError) as exc:
                    raise InterpolatorFormatError(
                        "malformed", f"malformed interpolator payload: {exc}"
                    ) from exc
        elif data is not None:
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
            if saved and interpolator._frame_at_data is None:
                interpolator._note_frame(1)

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
        element_volume: float | None = None,
        support=None,
        value_constraints: np.ndarray | None = None,
        gradient_norm_constraints: np.ndarray | None = None,
        gradient_constraints: np.ndarray | None = None,
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
