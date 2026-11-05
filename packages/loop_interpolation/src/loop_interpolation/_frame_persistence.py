"""Persistence of the local coordinate frame of an interpolator.

This module stores the frame of an interpolator in a payload and reads it
back.
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np
from loop_common.geometry import BoundingBox, RigidRotationError, validate_rigid_rotation
from loop_common.geometry._bounding_box import LoopValueError
from loop_common.utils import LoopValueError as UtilsLoopValueError

_PARSE_ERRORS = (
    RigidRotationError,
    LoopValueError,
    UtilsLoopValueError,
    KeyError,
    TypeError,
    ValueError,
)

FORMAT_VERSION: int = 2
"""Version of the frame payload format."""

CONSTRAINT_SPACE: str = "local"
"""Space in which the constraints of a frame are stored."""

REASONS: tuple[str, ...] = (
    "legacy_no_transform",
    "unsupported_version",
    "malformed",
    "conflict",
)
"""Valid reasons for an `InterpolatorFormatError`."""


class InterpolatorFormatError(ValueError):
    """Error raised when a stored interpolator payload cannot be used.

    Parameters
    ----------
    reason : str
        Short code for the cause. It must be one of `REASONS`.
    message : str
        Text that describes the error.

    Raises
    ------
    ValueError
        If `reason` is not one of `REASONS`.
    """

    def __init__(self, reason: str, message: str) -> None:
        if reason not in REASONS:
            raise ValueError(f"reason must be one of {REASONS}, got {reason!r}")
        super().__init__(message)
        self.reason = reason


class StaleFrameError(ValueError):
    """Error raised when the local transform changed after constraints were stored."""


def frame_payload(bounding_box: BoundingBox | None, dimensions: int) -> dict:
    """Build the payload that stores a frame.

    Parameters
    ----------
    bounding_box : BoundingBox or None
        Bounding box of the frame. If None, the frame is the identity.
    dimensions : int
        Number of spatial dimensions.

    Returns
    -------
    dict
        Payload with the keys `format_version`, `constraint_space`,
        `transform` and `bounding_box`, in that order.
    """
    if bounding_box is None:
        local_origin = np.zeros(dimensions)
        local_rotation = np.eye(dimensions)
        box = None
    else:
        local_origin = bounding_box.local_origin
        local_rotation = bounding_box.local_rotation
        box = {
            key: value
            for key, value in bounding_box.to_dict().items()
            if key not in ("local_origin", "local_rotation")
        }
    return {
        "format_version": FORMAT_VERSION,
        "constraint_space": CONSTRAINT_SPACE,
        "transform": {
            "dimensions": dimensions,
            "local_origin": local_origin.tolist(),
            "local_rotation": local_rotation.tolist(),
        },
        "bounding_box": box,
    }


def payload_kind(payload: dict, caller_frame: bool) -> str:
    """Return the kind of a stored payload.

    Parameters
    ----------
    payload : dict
        Stored payload.
    caller_frame : bool
        True if the caller gives a frame of its own.

    Returns
    -------
    str
        `"v2"` or `"legacy"`.

    Raises
    ------
    InterpolatorFormatError
        If the payload has no transform and the caller gives no frame
        (`legacy_no_transform`), if the version is above `FORMAT_VERSION`
        (`unsupported_version`), or if the payload is not valid (`malformed`).
    """
    version = payload.get("format_version")
    if version is None:
        if any(key in payload for key in ("constraint_space", "transform", "bounding_box")):
            raise InterpolatorFormatError(
                "malformed",
                "malformed interpolator payload: legacy payload has frame keys "
                "but no format_version",
            )
        if not caller_frame:
            raise InterpolatorFormatError(
                "legacy_no_transform",
                "legacy interpolator payload has no transform: pass local_origin "
                "and local_rotation, or declare_identity=True",
            )
        return "legacy"
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        raise InterpolatorFormatError(
            "malformed",
            f"malformed interpolator payload: format_version {version!r} is not "
            "an integer of 1 or more",
        )
    if version > FORMAT_VERSION:
        raise InterpolatorFormatError(
            "unsupported_version",
            f"unsupported interpolator format_version {version}: this reader "
            f"supports version {FORMAT_VERSION}",
        )
    if version == 1:
        raise InterpolatorFormatError(
            "malformed",
            "malformed interpolator payload: no reader writes format_version 1",
        )
    return "v2"


def _malformed(detail: str) -> InterpolatorFormatError:
    """Make a `malformed` error with the standard prefix."""
    return InterpolatorFormatError("malformed", f"malformed interpolator payload: {detail}")


def _parse_transform(
    origin: Any, rotation: Any, dimensions: int, name: str
) -> tuple[np.ndarray, np.ndarray]:
    """Check an origin and a rotation and return them as float arrays.

    Parameters
    ----------
    origin : Any
        Candidate origin.
    rotation : Any
        Candidate rotation matrix.
    dimensions : int
        Number of spatial dimensions.
    name : str
        Text that tells where the values come from, for messages.

    Returns
    -------
    tuple of np.ndarray
        The origin and the rotation as float arrays.

    Raises
    ------
    InterpolatorFormatError
        With reason `malformed`, if a value is not valid.
    """
    try:
        origin_array = np.array(origin, dtype=float)
        if origin_array.shape != (dimensions,):
            raise _malformed(
                f"{name} local_origin must have shape ({dimensions},), got {origin_array.shape}"
            )
        if not np.all(np.isfinite(origin_array)):
            raise _malformed(f"{name} local_origin must be finite, got {origin_array.tolist()}")
        rotation_array = validate_rigid_rotation(rotation, dimensions)
    except InterpolatorFormatError:
        raise
    except _PARSE_ERRORS as exc:
        raise _malformed(f"{name} transform is not valid: {exc}") from exc
    return origin_array, rotation_array


def resolve_frame(
    payload: dict,
    kind: str,
    support: object,
    dimensions: int,
    *,
    local_origin: Any = None,
    local_rotation: Any = None,
    declare_identity: bool = False,
) -> BoundingBox | None:
    """Resolve the frame of a stored payload.

    The function reads `payload`. It does not change `payload` or any object
    of the caller.

    Parameters
    ----------
    payload : dict
        Frame keys of the stored payload: `constraint_space`, `transform`
        and `bounding_box`. It can also hold `format_version`.
    kind : str
        Kind of the payload, from `payload_kind`: `"v2"` or `"legacy"`.
    support : object
        Support of the interpolator. It must have a `nodes` array.
    dimensions : int
        Number of spatial dimensions.
    local_origin : array_like, optional
        World-space origin of the local frame, given by the caller.
    local_rotation : array_like, optional
        Rotation matrix of the local frame, given by the caller.
    declare_identity : bool, optional
        True if the caller declares an identity frame.

    Returns
    -------
    BoundingBox or None
        Bounding box to attach. None for an identity frame with no box.

    Raises
    ------
    InterpolatorFormatError
        With reason `conflict` if the inputs disagree with each other or with
        the saved transform. With reason `malformed` if a value is not valid.
    """
    # Step 1: declare_identity cannot go with a transform value.
    if declare_identity and (local_origin is not None or local_rotation is not None):
        given = [
            name
            for name, value in (("local_origin", local_origin), ("local_rotation", local_rotation))
            if value is not None
        ]
        raise InterpolatorFormatError(
            "conflict",
            "conflicting transform: declare_identity=True cannot be used with "
            + " and ".join(given),
        )
    caller_frame = local_origin is not None or local_rotation is not None or declare_identity

    # Step 2: check the transform of the caller.
    origin = rotation = None
    if caller_frame:
        origin_value = (
            np.zeros(dimensions) if declare_identity or local_origin is None else local_origin
        )
        rotation_value = (
            np.eye(dimensions) if declare_identity or local_rotation is None else local_rotation
        )
        origin, rotation = _parse_transform(origin_value, rotation_value, dimensions, "supplied")

    # Step 3: a legacy payload has no saved frame.
    if kind == "legacy":
        if not caller_frame:
            raise _malformed("legacy payload has no transform and the caller gives none")
        if declare_identity:
            return None
        try:
            nodes = np.asarray(support.nodes, dtype=float)
        except (AttributeError, TypeError, ValueError) as exc:
            raise _malformed(f"support has no usable nodes array: {exc}") from exc
        if nodes.ndim != 2 or nodes.shape[1] != dimensions or nodes.shape[0] == 0:
            raise _malformed(f"support.nodes must have shape (n, {dimensions}), got {nodes.shape}")
        return bounding_box_from_support(support, dimensions, origin, rotation)
    if kind != "v2":
        raise _malformed(f"unknown payload kind {kind!r}")

    # Step 4: read the saved frame of a v2 payload.
    space = payload.get("constraint_space")
    if space != CONSTRAINT_SPACE:
        raise _malformed(f"constraint_space must be {CONSTRAINT_SPACE!r}, got {space!r}")
    transform = payload.get("transform")
    if not isinstance(transform, dict):
        raise _malformed(f"transform must be a dict, got {transform!r}")
    try:
        saved_dimensions = transform["dimensions"]
        saved_origin_value = transform["local_origin"]
        saved_rotation_value = transform["local_rotation"]
    except KeyError as exc:
        raise _malformed(f"transform has no key {exc}") from exc
    if (
        not isinstance(saved_dimensions, int)
        or isinstance(saved_dimensions, bool)
        or saved_dimensions != dimensions
    ):
        raise _malformed(
            f"transform dimensions {saved_dimensions!r} does not equal the interpolator "
            f"dimensions {dimensions}"
        )
    saved_origin, saved_rotation = _parse_transform(
        saved_origin_value, saved_rotation_value, dimensions, "saved"
    )
    box_data = payload.get("bounding_box")
    if box_data is not None and not isinstance(box_data, dict):
        raise _malformed(f"bounding_box must be a dict or None, got {type(box_data).__name__}")
    if box_data is None and not (
        np.array_equal(saved_origin, np.zeros(dimensions))
        and np.array_equal(saved_rotation, np.eye(dimensions))
    ):
        raise _malformed("bounding_box is None but the saved transform is not the identity")

    # Step 5: the transform of the caller must equal the saved one.
    if caller_frame:
        same_origin = bool(np.allclose(origin, saved_origin, rtol=0, atol=1e-12))
        same_rotation = bool(np.allclose(rotation, saved_rotation, rtol=0, atol=1e-12))
        if not (same_origin and same_rotation):
            if declare_identity:
                detail = (
                    "declare_identity=True conflicts with the saved transform "
                    f"(saved local_origin={saved_origin.tolist()}, "
                    f"saved local_rotation={saved_rotation.tolist()})"
                )
            else:
                parts = []
                if not same_origin:
                    parts.append(
                        f"local_origin: saved {saved_origin.tolist()}, supplied {origin.tolist()}"
                    )
                if not same_rotation:
                    parts.append(
                        f"local_rotation: saved {saved_rotation.tolist()}, "
                        f"supplied {rotation.tolist()}"
                    )
                detail = "; ".join(parts)
            raise InterpolatorFormatError("conflict", f"conflicting transform: {detail}")

    # Step 6: build the saved box.
    if box_data is None:
        return None
    try:
        box = BoundingBox.from_dict(
            {
                **box_data,
                "local_origin": saved_origin.tolist(),
                "local_rotation": saved_rotation.tolist(),
            }
        )
    except _PARSE_ERRORS as exc:
        raise _malformed(f"bounding_box is not valid: {exc}") from exc
    if box.dimensions != dimensions:
        raise _malformed(f"bounding_box has {box.dimensions} dimensions, expected {dimensions}")
    return box


def bounding_box_from_support(
    support: object,
    dimensions: int,
    local_origin: np.ndarray,
    local_rotation: np.ndarray,
) -> BoundingBox:
    """Compute the bounding box of a support in the world frame.

    The box is an envelope of the support. Only `local_enclosure` and display
    use it. Evaluation never uses it.

    Parameters
    ----------
    support : object
        Support of the interpolator. It must have a `nodes` array in local
        coordinates.
    dimensions : int
        Number of spatial dimensions.
    local_origin : np.ndarray
        World-space origin of the local frame.
    local_rotation : np.ndarray
        Rotation matrix of the local frame.

    Returns
    -------
    BoundingBox
        World-space box with the local transform set.
    """
    nodes = np.asarray(support.nodes, dtype=float)
    local_min = nodes.min(axis=0)
    local_max = nodes.max(axis=0)
    corners = np.array(list(itertools.product(*zip(local_min, local_max))))
    origin = np.asarray(local_origin, dtype=float)
    rotation = np.asarray(local_rotation, dtype=float)
    world = corners @ rotation + origin
    box = BoundingBox(dimensions=dimensions, origin=world.min(axis=0), maximum=world.max(axis=0))
    box.set_local_transform(local_origin=origin, rotation_matrix=rotation)
    return box
