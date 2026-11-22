"""Filter and colour commands: input checks, colour settings, results and requests."""

from __future__ import annotations

import math
import numbers
import threading
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from loop_link import encode
from loop_link.errors import InvalidInputError
from loop_link.protocol import Frame, ProtocolError

if TYPE_CHECKING:
    from loop_link.tool import ToolConnection

__all__ = [
    "ColourSettings",
    "FilterResult",
    "check_codes",
    "check_number",
    "check_rgb",
    "check_text",
    "format_number",
    "list_colour_map_names",
    "parse_colour_map_names",
    "raise_for_command_result",
    "run_command",
]

_FEATURE = "filter and colour commands (geoviewer protocol 1.11)"
_SCALES = ("linear", "log")
_CODE_MAX = 2**32 - 1


def format_number(value: numbers.Real) -> str:
    """Return the text of a number, as the ``repr`` of its float value.

    Parameters
    ----------
    value : number
        The number to format.

    Returns
    -------
    str
        The text of the float value.
    """
    return repr(float(value))


def check_number(value: object, name: str) -> float:
    """Check that a value is a finite number and return it as a float.

    Parameters
    ----------
    value : number
        The value to check. ``bool``, NaN and infinity are not valid.
    name : str
        The name of the argument, for the error message.

    Returns
    -------
    float
        The value as a float.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The value is not a finite number.
    """
    if (
        isinstance(value, (bool, str))
        or not isinstance(value, numbers.Real)
        or not math.isfinite(float(value))
    ):
        raise InvalidInputError(f"'{name}' must be a finite number. Got {value!r}.")
    return float(value)


def check_text(value: object, name: str) -> str:
    """Check that a value is a text that is not empty.

    Parameters
    ----------
    value : str
        The value to check.
    name : str
        The name of the argument, for the error message.

    Returns
    -------
    str
        The text.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The value is not a ``str`` or it is empty.
    """
    if not isinstance(value, str) or not value:
        raise InvalidInputError(f"'{name}' must be a text that is not empty. Got {value!r}.")
    return value


def _check_int(value: object, name: str, low: int, high: int) -> int:
    """Return a value as an ``int`` when it is an integer in a range."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise InvalidInputError(f"'{name}' must contain only integers. Got {value!r}.")
    number = int(value)
    if not low <= number <= high:
        raise InvalidInputError(f"'{name}' must contain only integers from {low} to {high}.")
    return number


def check_rgb(value: object, name: str) -> tuple[int, int, int]:
    """Check that a value is three integers from 0 to 255.

    Parameters
    ----------
    value : sequence of int
        The value to check.
    name : str
        The name of the argument, for the error message.

    Returns
    -------
    tuple of int
        The red, green and blue values.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The value is not three integers from 0 to 255.
    """
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise InvalidInputError(f"'{name}' must be three integers from 0 to 255. Got {value!r}.")
    red, green, blue = (_check_int(item, name, 0, 255) for item in value)
    return (red, green, blue)


def check_codes(values: Iterable[Any], name: str) -> list[int]:
    """Check that each value is an integer code from 0 to 2**32-1.

    Parameters
    ----------
    values : iterable of int
        The codes to check. A numpy array is valid.
    name : str
        The name of the argument, for the error message.

    Returns
    -------
    list of int
        The codes as Python integers.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The value is not iterable or one code is not valid.
    """
    if isinstance(values, (str, bytes)):
        raise InvalidInputError(f"'{name}' must be a list of integer codes. Got {values!r}.")
    try:
        items = list(values)
    except TypeError as error:
        raise InvalidInputError(f"'{name}' must be a list of integer codes.") from error
    return [_check_int(item, name, 0, _CODE_MAX) for item in items]


@dataclass(frozen=True)
class ColourSettings:
    """The colour settings of one object. A value of ``None`` means no change.

    Parameters
    ----------
    array : str, optional
        The name of the array to colour by.
    colour_map : str, optional
        The name of the colour map.
    range : tuple of float, optional
        The low and high limit of the colour range. The low limit is less than the high limit.
    scale : str, optional
        The colour scale, ``"linear"`` or ``"log"``.
    plain : tuple of int, optional
        The plain colour as red, green and blue values from 0 to 255.

    Raises
    ------
    loop_link.errors.InvalidInputError
        A value is not valid. The message names the field.
    """

    array: str | None = None
    colour_map: str | None = None
    range: tuple[float, float] | None = None
    scale: str | None = None
    plain: tuple[int, int, int] | None = None

    def __post_init__(self) -> None:
        """Check the values and store them in their final form."""
        if self.array is not None:
            check_text(self.array, "array")
        if self.colour_map is not None:
            check_text(self.colour_map, "colour_map")
        if self.range is not None:
            limits = self.range
            if isinstance(limits, (str, bytes)) or not isinstance(limits, Sequence):
                raise InvalidInputError(f"'range' must be two numbers. Got {limits!r}.")
            if len(limits) != 2:
                raise InvalidInputError(f"'range' must be two numbers. Got {limits!r}.")
            low = check_number(limits[0], "range")
            high = check_number(limits[1], "range")
            if not low < high:
                raise InvalidInputError("'range' must have a low limit less than the high limit.")
            object.__setattr__(self, "range", (low, high))
        if self.scale is not None and self.scale not in _SCALES:
            raise InvalidInputError(f"'scale' must be 'linear' or 'log'. Got {self.scale!r}.")
        if self.plain is not None:
            object.__setattr__(self, "plain", check_rgb(self.plain, "plain"))


def _count(result: Mapping[str, Any], key: str) -> int | None:
    """Return a count of a result, or ``None`` when it is absent."""
    value = result.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ProtocolError(f"The viewer sent a filter result '{key}' that is not valid.")
    return value


@dataclass(frozen=True)
class FilterResult:
    """The result of a filter call. A field is ``None`` when the viewer gives no value.

    Parameters
    ----------
    passing : int, optional
        The number of elements that pass the filters.
    total : int, optional
        The number of elements of the object.
    state : str, optional
        ``"unbuilt"`` when the viewer did not build the filters, or ``"building"`` when
        the viewer stopped to wait and the build goes on.
    """

    passing: int | None = None
    total: int | None = None
    state: str | None = None

    @classmethod
    def from_result(cls, result: Mapping[str, Any]) -> FilterResult:
        """Read a ``command.result`` body of ``filter.edit``.

        Parameters
        ----------
        result : mapping
            The ``result`` body, with the status ``done``.

        Returns
        -------
        FilterResult
            The counts and the state of the body.

        Raises
        ------
        loop_link.protocol.ProtocolError
            A count is not a whole number of 0 or more, or the state is not a text.
        """
        state = result.get("state")
        if state is not None and not isinstance(state, str):
            raise ProtocolError("The viewer sent a filter result 'state' that is not a text.")
        return cls(_count(result, "passing"), _count(result, "total"), state)


def raise_for_command_result(result: Mapping[str, Any], *, command: str) -> None:
    """Raise an error when a ``command.result`` body is not ``done``.

    Parameters
    ----------
    result : mapping
        The ``result`` body of the frame.
    command : str
        The name of the command, for the error message.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The status is ``refused``.
    loop_link.protocol.ProtocolError
        The status is not known.
    """
    status = result.get("status")
    if status == "done":
        return
    if status == "refused":
        raise InvalidInputError(
            f"The viewer refused the command '{command}': {result.get('text', '')}"
        )
    raise ProtocolError(f"The viewer sent an unknown result status '{status}'.")


def parse_colour_map_names(frame: Frame) -> tuple[str, ...] | None:
    """Read the colour map names from a ``commands`` frame.

    Parameters
    ----------
    frame : loop_link.protocol.Frame
        The ``commands`` frame of the viewer.

    Returns
    -------
    tuple of str or None
        The names, or ``None`` when the frame has no ``colour.set_map`` entry or the entry
        has no list of names.

    Raises
    ------
    loop_link.protocol.ProtocolError
        The list of names is not a list of ``str``.
    """
    entries = frame.header.get("commands")
    if not isinstance(entries, list):
        return None
    for entry in entries:
        if not isinstance(entry, Mapping) or entry.get("name") != "colour.set_map":
            continue
        schema = entry.get("input_schema")
        properties = schema.get("properties") if isinstance(schema, Mapping) else None
        name = properties.get("name") if isinstance(properties, Mapping) else None
        if not isinstance(name, Mapping) or "enum" not in name:
            return None
        names = name["enum"]
        if not isinstance(names, list) or not all(isinstance(item, str) for item in names):
            raise ProtocolError("The viewer sent colour map names that are not a list of text.")
        return tuple(names)
    return None


def run_command(
    tool: ToolConnection,
    command: str,
    input: Mapping[str, Any],
    *,
    object_id: str | None = None,
    cancel: threading.Event | None = None,
    timeout: float | None = None,
    progress: Callable[[float], None] | None = None,
) -> Mapping[str, Any]:
    """Send one command and return its result.

    Parameters
    ----------
    tool : loop_link.tool.ToolConnection
        The tool connection.
    command : str
        The name of the command.
    input : mapping
        The input of the command.
    object_id : str, optional
        The identifier of the object that the command is about, for error messages.
    cancel : threading.Event, optional
        When it is set, the request is cancelled.
    timeout : float, optional
        The time in seconds to wait for the next frame.
    progress : callable, optional
        Called with the fraction of each ``command.progress`` frame.

    Returns
    -------
    mapping
        The ``result`` body, with the status ``done``.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The viewer refused the command.
    loop_link.protocol.ProtocolError
        The result is not a mapping with a ``str`` status, or the status is not known.
    """
    frame = tool.request(
        lambda number: encode.command_frame(number, command, input),
        until=("command.result",),
        feature=_FEATURE,
        object_id=object_id,
        progress=progress,
        cancel=cancel,
        timeout=timeout,
    )
    result = frame.header.get("result")
    if not isinstance(result, Mapping) or not isinstance(result.get("status"), str):
        raise ProtocolError("The viewer sent a command result that is not valid.")
    raise_for_command_result(result, command=command)
    return result


def list_colour_map_names(tool: ToolConnection) -> tuple[str, ...] | None:
    """Ask the viewer for the names of its colour maps.

    Parameters
    ----------
    tool : loop_link.tool.ToolConnection
        The tool connection.

    Returns
    -------
    tuple of str or None
        The names, or ``None`` when the viewer gives no list.
    """
    frame = tool.request(
        lambda number: encode.commands_list_frame(number),
        until=("commands",),
        feature=_FEATURE,
    )
    return parse_colour_map_names(frame)
