"""Persistence of the local coordinate frame of an interpolator.

This module stores the frame of an interpolator in a payload and reads it
back.
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np
from loop_common.geometry import BoundingBox

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


def resolve_frame(
    payload: dict,
    kind: str,
    support: object,
    dimensions: int,
    *,
    local_origin: Any = None,
    local_rotation: Any = None,
    declare_identity: bool = False,
) -> Any:
    """Resolve the frame of a stored payload.

    Parameters
    ----------
    payload : dict
        Stored payload.
    kind : str
        Kind of the payload, from `payload_kind`.
    support : object
        Support of the interpolator.
    dimensions : int
        Number of spatial dimensions.
    local_origin : Any, optional
        Origin of the local frame.
    local_rotation : Any, optional
        Rotation of the local frame.
    declare_identity : bool, optional
        True if the caller declares an identity frame.

    Returns
    -------
    Any
        Resolved frame.
    """
    raise NotImplementedError("S3 stub")


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
