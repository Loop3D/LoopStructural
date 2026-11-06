from __future__ import annotations

import json
import os
import struct
from pathlib import Path

import numpy as np
import pytest

from loop_link import protocol
from loop_link.protocol import (
    MAX_FRAME_BYTES,
    MAX_HEADER_BYTES,
    PROTOCOL_VERSION,
    Frame,
    ProtocolError,
    decode_frame,
    decode_header,
    encode_frame,
)

DTYPES = {
    "f32": np.float32,
    "f64": np.float64,
    "u32": np.uint32,
    "i32": np.int32,
    "u8": np.uint8,
}


def test_constants():
    # The limits and the version are fixed by the protocol.
    assert MAX_FRAME_BYTES == 256 * 1024 * 1024
    assert MAX_HEADER_BYTES == 1024 * 1024
    assert PROTOCOL_VERSION == {"major": 1, "minor": 14}


@pytest.mark.parametrize("name", list(DTYPES))
def test_round_trip_each_dtype(name):
    # Encode one array of each dtype, then decode it. The data must stay the same.
    array = np.arange(6).reshape(2, 3).astype(DTYPES[name])
    data = encode_frame({"type": "test"}, [("a", array)])
    frame = decode_frame(data)
    assert isinstance(frame, Frame)
    assert frame.header["type"] == "test"
    assert frame.header["version"] == PROTOCOL_VERSION
    assert frame.header["buffers"] == [{"name": "a", "dtype": name, "shape": [2, 3]}]
    assert frame.buffers["a"].dtype == array.dtype
    assert frame.buffers["a"].shape == (2, 3)
    np.testing.assert_array_equal(frame.buffers["a"], array)


def test_round_trip_several_buffers():
    # A frame can hold more than one buffer. Each buffer keeps its name.
    a = np.arange(4, dtype=np.float64)
    b = np.arange(3, dtype=np.uint8)
    frame = decode_frame(encode_frame({"type": "test"}, [("a", a), ("b", b)]))
    assert set(frame.buffers) == {"a", "b"}
    np.testing.assert_array_equal(frame.buffers["a"], a)
    np.testing.assert_array_equal(frame.buffers["b"], b)


def test_zero_element_array_is_valid():
    # An array with no elements is a valid buffer.
    array = np.zeros((0, 3), dtype=np.float32)
    frame = decode_frame(encode_frame({"type": "test"}, [("empty", array)]))
    assert frame.buffers["empty"].shape == (0, 3)
    assert frame.buffers["empty"].dtype == np.float32


def test_no_buffers():
    # A frame with no buffers has an empty buffer list.
    frame = decode_frame(encode_frame({"type": "heartbeat", "seq": 1}))
    assert frame.header["buffers"] == []
    assert frame.buffers == {}


def test_decode_header_only():
    # decode_header reads only the header.
    data = encode_frame({"type": "test"}, [("a", np.zeros(2, dtype=np.float32))])
    header = decode_header(data)
    assert header["type"] == "test"
    assert header["version"] == PROTOCOL_VERSION


def test_wrong_frame_length_is_error():
    # A frame with extra or missing bytes does not match 4 + header_len + sum(sizes).
    data = encode_frame({"type": "test"}, [("a", np.zeros(4, dtype=np.float32))])
    with pytest.raises(ProtocolError):
        decode_frame(data + b"\x00")
    with pytest.raises(ProtocolError):
        decode_frame(data[:-1])


def test_oversize_header_is_error_on_encode():
    # A header over 1 MiB is not allowed.
    with pytest.raises(ProtocolError):
        encode_frame({"type": "test", "pad": "x" * (MAX_HEADER_BYTES + 1)})


@pytest.mark.parametrize("size", [0, 1, 3])
def test_frame_under_four_bytes_is_error(size):
    # A frame needs 4 bytes for the header length.
    with pytest.raises(ProtocolError):
        decode_frame(b"\x00" * size)


def test_nan_in_header_is_error():
    # JSON has no NaN. The encoder must refuse it.
    with pytest.raises((ProtocolError, ValueError)):
        encode_frame({"type": "test", "value": float("nan")})


def test_unsupported_dtype_is_error():
    # Only f32, f64, u32, i32 and u8 are allowed.
    array = np.zeros((2, 3), dtype=np.complex128)
    with pytest.raises(ProtocolError):
        encode_frame({"type": "test"}, [("a", array)])


def test_byte_literal_layout():
    # The frame starts with the header length (u32, little-endian).
    # The JSON header has no spaces.
    data = encode_frame({"type": "heartbeat", "seq": 1})
    json_part = data[4:]
    assert data[:4] == struct.pack("<I", len(json_part))
    assert b" " not in json_part
    header = json.loads(json_part)
    assert isinstance(header, dict)
    assert header["type"] == "heartbeat"
    assert header["seq"] == 1


def test_frame_over_max_is_error_on_encode(monkeypatch):
    # A frame over MAX_FRAME_BYTES is not allowed.
    monkeypatch.setattr(protocol, "MAX_FRAME_BYTES", 100)
    with pytest.raises(ProtocolError):
        encode_frame({"type": "test"}, [("a", np.zeros(100, dtype=np.float64))])


def test_big_endian_input_gives_little_endian_bytes():
    # An input in big-endian order must be sent in little-endian order.
    array = np.array([1.0, 2.0], dtype=">f4")
    frame = decode_frame(encode_frame({"type": "test"}, [("a", array)]))
    assert frame.buffers["a"].tobytes() == np.array([1.0, 2.0], "<f4").tobytes()


def test_non_contiguous_input_gives_contiguous_bytes():
    # A non-contiguous input must be sent as contiguous bytes.
    array = np.array([1.0, 0.0, 2.0, 0.0], dtype="<f4")[::2]
    assert not array.flags["C_CONTIGUOUS"]
    frame = decode_frame(encode_frame({"type": "test"}, [("a", array)]))
    assert frame.buffers["a"].tobytes() == np.array([1.0, 2.0], "<f4").tobytes()


def test_header_length_over_data_is_error():
    # The header length field is larger than the data.
    data = struct.pack("<I", 1000) + b"{}"
    with pytest.raises(ProtocolError):
        decode_frame(data)
    with pytest.raises(ProtocolError):
        decode_header(data)


def _raw_frame(buffers, payload=b""):
    # Build a frame with a hand-written buffer list.
    body = json.dumps({"type": "test", "buffers": buffers}).encode("utf-8")
    return struct.pack("<I", len(body)) + body + payload


@pytest.mark.parametrize(
    "entry",
    [
        {"name": "a", "dtype": "u8", "shape": [4294967296, 4294967296]},
        {"name": "a", "dtype": "u8", "shape": [3037000500, 3037000500]},
        {"name": ["x"], "dtype": "u8", "shape": [1]},
        {"name": "a", "dtype": "u8", "shape": [1.0]},
        {"name": "a", "dtype": "u8", "shape": [True]},
        {"name": "a", "dtype": "u8", "shape": [-1]},
        {"name": "a", "dtype": "nope", "shape": [1]},
    ],
)
def test_bad_buffer_description_is_protocol_error(entry):
    # A bad header must raise ProtocolError and no other error.
    with pytest.raises(ProtocolError):
        decode_frame(_raw_frame([entry], b"\x00"))


def test_duplicate_buffer_name_is_error():
    # Two buffers with the same name are not allowed.
    entry = {"name": "a", "dtype": "u8", "shape": [1]}
    with pytest.raises(ProtocolError):
        decode_frame(_raw_frame([entry, entry], b"\x00\x00"))


def test_frame_over_max_is_error_on_decode(monkeypatch):
    # A frame over MAX_FRAME_BYTES is not allowed.
    data = encode_frame({"type": "test"})
    monkeypatch.setattr(protocol, "MAX_FRAME_BYTES", len(data) - 1)
    with pytest.raises(ProtocolError):
        decode_frame(data)


_FIXTURE_DIR = os.environ.get("LOOP_LINK_GEOVIEWER_FIXTURES")


@pytest.mark.skipif(not _FIXTURE_DIR, reason="LOOP_LINK_GEOVIEWER_FIXTURES is not set")
def test_geoviewer_fixtures_decode():
    # Each fixture file from the viewer must decode without error.
    files = sorted(Path(_FIXTURE_DIR).glob("*.bin"))
    assert files, "no *.bin files in the fixture directory"
    for path in files:
        decode_frame(path.read_bytes())
