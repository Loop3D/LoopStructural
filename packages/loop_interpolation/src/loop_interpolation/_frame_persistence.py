"""Persistence of the local coordinate frame of an interpolator.

This module stores the frame of an interpolator in a payload and reads it
back. The stubs in this module are filled in by a later step of the spec.
"""

from __future__ import annotations

from typing import Any

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
    """Error raised when a stored frame no longer matches its support."""


def frame_payload(bounding_box: Any, dimensions: int) -> dict:
    """Build the payload that stores a frame.

    Parameters
    ----------
    bounding_box : Any
        Bounding box of the frame.
    dimensions : int
        Number of spatial dimensions.

    Returns
    -------
    dict
        Payload that can be stored.
    """
    raise NotImplementedError("S3 stub")


def payload_kind(payload: dict, caller_frame: bool) -> str:
    """Return the kind of a stored payload.

    Parameters
    ----------
    payload : dict
        Stored payload.
    caller_frame : bool
        True if the caller asks for a frame of its own.

    Returns
    -------
    str
        Kind of the payload.
    """
    raise NotImplementedError("S3 stub")


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
    local_origin: Any,
    local_rotation: Any,
) -> Any:
    """Compute the bounding box of a support in the local frame.

    Parameters
    ----------
    support : object
        Support of the interpolator.
    dimensions : int
        Number of spatial dimensions.
    local_origin : Any
        Origin of the local frame.
    local_rotation : Any
        Rotation of the local frame.

    Returns
    -------
    Any
        Bounding box in the local frame.
    """
    raise NotImplementedError("S3 stub")
