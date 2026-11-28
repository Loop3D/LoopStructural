"""Error types that loop_link raises.

All errors derive from :class:`LoopLinkError`, so a caller can catch every
loop_link error with one ``except`` clause.
"""

from __future__ import annotations

from dataclasses import dataclass


class LoopLinkError(Exception):
    """Base class for all loop_link errors."""


class NoViewerError(LoopLinkError):
    """No viewer is available.

    Parameters
    ----------
    message : str
        The error text. The default text tells the user how to start a viewer.
    """

    def __init__(
        self,
        message: str = "No viewer is running. Start a viewer before you send objects.",
    ) -> None:
        super().__init__(message)


class TokenError(LoopLinkError):
    """The access token is missing, wrong or expired.

    Parameters
    ----------
    message : str
        The error text. The default text tells the user how to get a new token.
    """

    def __init__(
        self,
        message: str = "The access token is not valid. Create a new token and try again.",
    ) -> None:
        super().__init__(message)


class BusyError(LoopLinkError):
    """The viewer is busy with another request.

    Parameters
    ----------
    message : str
        The error text. The default text tells the user to wait and retry.
    """

    def __init__(
        self,
        message: str = "The viewer is busy. Wait a moment and try again.",
    ) -> None:
        super().__init__(message)


class InvalidInputError(LoopLinkError):
    """An input value is not valid.

    Parameters
    ----------
    message : str
        The error text. The default text tells the user to check the input.
    """

    def __init__(
        self,
        message: str = "An input value is not valid. Check the input and try again.",
    ) -> None:
        super().__init__(message)


class ViewerClosedError(LoopLinkError):
    """The viewer was closed.

    Parameters
    ----------
    message : str
        The error text. The default text tells the user to open the viewer again.
    """

    def __init__(
        self,
        message: str = "The viewer was closed. Open the viewer again and retry.",
    ) -> None:
        super().__init__(message)


class MissingFeatureError(LoopLinkError):
    """A required feature is not available.

    Parameters
    ----------
    feature : str
        The name of the missing feature.

    Attributes
    ----------
    feature : str
        The name of the missing feature.
    """

    def __init__(self, feature: str) -> None:
        self.feature = feature
        super().__init__(
            f"The feature '{feature}' is not available in this installation. "
            f"Check that the '{feature}' feature is installed and enabled, then retry."
        )


class ObjectTooLargeError(LoopLinkError):
    """An object is larger than the allowed limit.

    Parameters
    ----------
    object_id : str
        The identifier of the object.
    size : int
        The size of the object.
    limit : int
        The largest size that is allowed.

    Attributes
    ----------
    object_id : str
        The identifier of the object.
    size : int
        The size of the object.
    limit : int
        The largest size that is allowed.
    """

    def __init__(self, object_id: str, size: int, limit: int) -> None:
        self.object_id = object_id
        self.size = size
        self.limit = limit
        super().__init__(
            f"Object '{object_id}' is {size} units, which is larger than the limit of "
            f"{limit} units. Reduce the size of the object or split it into smaller "
            "objects, then try again."
        )


class StaleRevisionError(LoopLinkError):
    """An object is not at the revision that the caller expects.

    Parameters
    ----------
    object_id : str
        The identifier of the object.
    revision : int, optional
        The revision that the caller expected.

    Attributes
    ----------
    object_id : str
        The identifier of the object.
    revision : int, optional
        The revision that the caller expected.
    """

    def __init__(self, object_id: str, revision: int | None = None) -> None:
        self.object_id = object_id
        self.revision = revision
        super().__init__(
            f"The object '{object_id}' is not at that revision in the viewer. "
            "Read the selection or fetch the object again, then try again."
        )


class CancelledError(LoopLinkError):
    """The request was cancelled.

    Parameters
    ----------
    message : str
        The error text. The default text says that the request was cancelled.
    """

    def __init__(self, message: str = "The request was cancelled.") -> None:
        super().__init__(message)


class UnsupportedObjectError(LoopLinkError):
    """The object type cannot be fetched by loop_link.

    Parameters
    ----------
    object_id : str
        The identifier of the object.
    object_type : str
        The type of the object.

    Attributes
    ----------
    object_id : str
        The identifier of the object.
    object_type : str
        The type of the object.
    """

    def __init__(self, object_id: str, object_type: str) -> None:
        self.object_id = object_id
        self.object_type = object_type
        super().__init__(
            f"The object '{object_id}' has the type '{object_type}'. loop_link can fetch "
            "Points, Curve, Surface, BlockModel, TetraMesh and Drillhole objects only."
        )


@dataclass(frozen=True)
class Refusal:
    """A reason why the viewer refused one object.

    Parameters
    ----------
    object_id : str
        The identifier of the refused object.
    code : str
        A short code that names the reason.
    message : str
        A message that explains the reason.
    """

    object_id: str
    code: str
    message: str


class ViewerRefusedError(LoopLinkError):
    """The viewer refused one or more objects.

    Parameters
    ----------
    refusals : list of Refusal
        The refusals, one for each refused object.

    Attributes
    ----------
    refusals : list of Refusal
        The refusals, one for each refused object.
    """

    def __init__(self, refusals: list[Refusal]) -> None:
        self.refusals = refusals
        count = len(refusals)
        super().__init__(
            f"The viewer refused {count} object(s). Read the code and message of each "
            "refusal, fix those objects, and try again."
        )


class MissingExtraError(LoopLinkError):
    """An optional extra package is not installed.

    Parameters
    ----------
    extra : str
        The name of the extra.

    Attributes
    ----------
    extra : str
        The name of the extra.
    """

    def __init__(self, extra: str) -> None:
        self.extra = extra
        super().__init__(
            f"The optional '{extra}' part of loop-link is not installed. "
            f"Run: pip install loop-link[{extra}]"
        )
