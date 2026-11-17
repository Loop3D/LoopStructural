"""The binary frame codec of the loop_link protocol.

A message is one binary frame. All numbers are little-endian::

    [u32 header_len][JSON header, UTF-8][buffer 0][buffer 1]...

The header lists the buffers in order. This module uses only numpy and the
standard library.
"""

from __future__ import annotations

import json
import math
import struct
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from .errors import LoopLinkError

PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 14}
TOOL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 18}
CANCEL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}
# The version of the commands.list and command.run frames.
COMMAND_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}

# The largest frame that is allowed, in bytes.
MAX_FRAME_BYTES: int = 256 * 1024 * 1024
# The largest JSON header that is allowed, in bytes.
MAX_HEADER_BYTES: int = 1024 * 1024
# The largest object that is allowed, in bytes. The encoder in encode.py uses it.
MAX_OBJECT_BYTES: int = 2 * 1024**3
# Beyond this distance from the origin an f32 coordinate is not exact to 1 mm.
F32_LIMIT: float = 4000.0

_LENGTH = struct.Struct("<I")

_DTYPE_NAMES: dict[str, str] = {
    np.dtype("<f4").str: "f32",
    np.dtype("<f8").str: "f64",
    np.dtype("<u4").str: "u32",
    np.dtype("<i4").str: "i32",
    np.dtype("u1").str: "u8",
}
_NAME_DTYPES: dict[str, np.dtype] = {
    "f32": np.dtype("<f4"),
    "f64": np.dtype("<f8"),
    "u32": np.dtype("<u4"),
    "i32": np.dtype("<i4"),
    "u8": np.dtype("u1"),
}


class ProtocolError(LoopLinkError):
    """A frame cannot be written or read."""


@dataclass(frozen=True)
class Frame:
    """A decoded frame.

    Parameters
    ----------
    header : dict
        The JSON header of the frame.
    buffers : dict of str to numpy.ndarray
        The arrays of the frame, by name. Each array has the shape from the header.
    """

    header: dict[str, Any]
    buffers: dict[str, np.ndarray]


def encode_frame(
    header: dict[str, Any],
    buffers: Sequence[tuple[str, np.ndarray]] = (),
    *,
    version: dict[str, int] | None = None,
) -> bytes:
    """Return the frame for a header and its named buffers.

    The function adds ``version`` and ``buffers`` to the header.

    Parameters
    ----------
    header : dict
        The header fields. The function does not change this dict.
    buffers : sequence of (str, numpy.ndarray)
        The named arrays. The dtype must be f32, f64, u32, i32 or u8.
    version : dict, optional
        The protocol version for the header. The default is ``PROTOCOL_VERSION``.
        The function does not change this dict.

    Returns
    -------
    bytes
        The frame.

    Raises
    ------
    ProtocolError
        If a dtype is not supported, or the header or the frame is too large.
    """
    descriptions = []
    payload = []
    for name, array in buffers:
        array = np.ascontiguousarray(array)
        little = array.dtype.newbyteorder("<")
        dtype_name = _DTYPE_NAMES.get(little.str)
        if dtype_name is None:
            raise ProtocolError(
                f"Buffer {name!r} has the type {array.dtype}, which is not supported. "
                "Use f32, f64, u32, i32 or u8."
            )
        descriptions.append({"name": name, "dtype": dtype_name, "shape": list(array.shape)})
        payload.append(array.astype(little, copy=False).tobytes())
    full = dict(header)
    full["version"] = PROTOCOL_VERSION if version is None else dict(version)
    full["buffers"] = descriptions
    header_bytes = json.dumps(full, separators=(",", ":"), allow_nan=False).encode("utf-8")
    if len(header_bytes) > MAX_HEADER_BYTES:
        raise ProtocolError(
            f"The header is {len(header_bytes)} bytes. The limit is {MAX_HEADER_BYTES} bytes."
        )
    total = _LENGTH.size + len(header_bytes) + sum(len(part) for part in payload)
    if total > MAX_FRAME_BYTES:
        raise ProtocolError(
            f"The frame is {total} bytes. The limit is {MAX_FRAME_BYTES} bytes. "
            "Send smaller objects."
        )
    return b"".join([_LENGTH.pack(len(header_bytes)), header_bytes, *payload])


def _read_header(data: bytes) -> tuple[dict[str, Any], int]:
    """Return the header and the offset where the buffers start."""
    if len(data) < _LENGTH.size:
        raise ProtocolError("The frame is shorter than 4 bytes.")
    (length,) = _LENGTH.unpack_from(data, 0)
    if length > MAX_HEADER_BYTES:
        raise ProtocolError(
            f"The header length is {length} bytes. The limit is {MAX_HEADER_BYTES} bytes."
        )
    end = _LENGTH.size + length
    if end > len(data):
        raise ProtocolError("The frame stops inside the header.")
    try:
        header = json.loads(data[_LENGTH.size : end].decode("utf-8"))
    except (ValueError, RecursionError) as error:
        raise ProtocolError(f"The header is not valid JSON: {error}") from error
    if not isinstance(header, dict):
        raise ProtocolError("The header is not a JSON object.")
    return header, end


def decode_header(data: bytes) -> dict[str, Any]:
    """Return the JSON header of a frame. The buffers are ignored.

    Parameters
    ----------
    data : bytes
        The frame.

    Returns
    -------
    dict
        The header.

    Raises
    ------
    ProtocolError
        If the frame is too short or the header is not valid.
    """
    return _read_header(data)[0]


def decode_frame(data: bytes) -> Frame:
    """Return the header and the arrays of a frame.

    Parameters
    ----------
    data : bytes
        The frame.

    Returns
    -------
    Frame
        The header and the named arrays.

    Raises
    ------
    ProtocolError
        If the frame is not valid, or its length does not match the header.
    """
    if len(data) > MAX_FRAME_BYTES:
        raise ProtocolError(
            f"The frame is {len(data)} bytes. The limit is {MAX_FRAME_BYTES} bytes."
        )
    header, offset = _read_header(data)
    descriptions = header.get("buffers", [])
    if not isinstance(descriptions, list):
        raise ProtocolError("The header field 'buffers' is not a list.")
    buffers: dict[str, np.ndarray] = {}
    for entry in descriptions:
        if not isinstance(entry, dict):
            raise ProtocolError(f"A buffer description is not valid: {entry!r}.")
        name = entry.get("name")
        dtype_name = entry.get("dtype")
        raw_shape = entry.get("shape")
        if not isinstance(name, str):
            raise ProtocolError(f"A buffer name is not a string: {name!r}.")
        if name in buffers:
            raise ProtocolError(f"The buffer name {name!r} is used more than once.")
        if not isinstance(dtype_name, str) or dtype_name not in _NAME_DTYPES:
            raise ProtocolError(f"Buffer {name!r} has an unknown dtype: {dtype_name!r}.")
        if not isinstance(raw_shape, list) or any(type(n) is not int or n < 0 for n in raw_shape):
            raise ProtocolError(f"Buffer {name!r} has a shape that is not valid: {raw_shape!r}.")
        dtype = _NAME_DTYPES[dtype_name]
        shape = tuple(raw_shape)
        count = math.prod(shape)
        size = count * dtype.itemsize
        if offset + size > len(data):
            raise ProtocolError(f"The frame stops inside buffer {name!r}.")
        try:
            buffers[name] = np.frombuffer(data, dtype=dtype, count=count, offset=offset).reshape(
                shape
            )
        except ValueError as error:
            raise ProtocolError(f"Buffer {name!r} cannot be read: {error}") from error
        offset += size
    if offset != len(data):
        raise ProtocolError(f"The frame is {len(data)} bytes. The header describes {offset} bytes.")
    return Frame(header=header, buffers=buffers)
