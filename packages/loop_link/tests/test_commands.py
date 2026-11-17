"""Tests for the command frames of the loop_link protocol."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from loop_link.commands import (
    ColourSettings,
    FilterResult,
    check_codes,
    check_number,
    check_rgb,
    check_text,
    format_number,
    list_colour_map_names,
    parse_colour_map_names,
    raise_for_command_result,
    run_command,
)
from loop_link.encode import command_frame, commands_list_frame
from loop_link.errors import InvalidInputError
from loop_link.protocol import Frame, ProtocolError, decode_frame

FEATURE = "filter and colour commands (geoviewer protocol 1.11)"


def test_commands_list_header():
    # The commands.list frame has the fixed header from the geoviewer fixtures.
    frame = decode_frame(commands_list_frame(1))
    assert frame.header == {
        "type": "commands.list",
        "version": {"major": 1, "minor": 11},
        "request": 1,
        "buffers": [],
    }


def test_command_frame_header():
    # The command.run frame carries the command name and the input.
    frame = decode_frame(command_frame(2, "object.set_visible", {"id": "p1", "visible": False}))
    assert frame.header == {
        "type": "command.run",
        "version": {"major": 1, "minor": 11},
        "request": 2,
        "command": "object.set_visible",
        "input": {"id": "p1", "visible": False},
        "buffers": [],
    }


def test_command_frames_have_no_buffers():
    # Neither frame carries binary buffers.
    assert decode_frame(commands_list_frame(1)).buffers == {}
    frame = decode_frame(command_frame(2, "object.set_visible", {"id": "p1", "visible": False}))
    assert frame.buffers == {}


def test_command_frame_does_not_change_input():
    # The caller's input mapping must keep its values after the call.
    input_ = {"id": "p1", "visible": False}
    before = dict(input_)
    command_frame(2, "object.set_visible", input_)
    assert input_ == before


def test_nested_input_round_trip():
    # A nested input comes back the same after encode and decode.
    input_ = {"object": "g", "edit": {"kind": "set_values", "array": "a", "values": [0, 2]}}
    frame = decode_frame(command_frame(3, "object.edit", input_))
    assert frame.header["input"] == input_


def test_nan_in_input_is_protocol_error():
    # JSON has no NaN. The encoder must refuse it.
    with pytest.raises(ProtocolError):
        command_frame(4, "object.set_visible", {"id": "p1", "value": float("nan")})


# The commands module: checks, settings, results and requests.


class StubTool:
    """A tool stub. It records the arguments of request and returns one frame."""

    def __init__(self, reply):
        self.reply = reply
        self.calls = []

    def request(self, build, **kwargs):
        self.calls.append((build, kwargs))
        return self.reply


def _result_frame(result):
    return Frame({"type": "command.result", "request": 1, "result": result}, {})


def _commands_frame(entries):
    return Frame({"type": "commands", "request": 1, "commands": entries}, {})


def _map_entry(enum):
    return {
        "name": "colour.set_map",
        "input_schema": {"properties": {"name": {"type": "string", "enum": enum}}},
    }


@pytest.mark.parametrize(
    ("value", "text"),
    [(0.5, "0.5"), (1e-05, "1e-05"), (3, "3.0"), (np.float32(0.5), "0.5")],
)
def test_format_number(value, text):
    # The text is the repr of the float value.
    assert format_number(value) == text


@pytest.mark.parametrize("value", [1, 2.5, np.int64(3), np.float64(1.5), np.float32(0.5)])
def test_check_number_accepts_numbers(value):
    # int, float and numpy numbers give a float.
    result = check_number(value, "x")
    assert isinstance(result, float)
    assert result == float(value)


@pytest.mark.parametrize(
    "value", [True, False, float("nan"), float("inf"), float("-inf"), "1", None]
)
def test_check_number_refuses(value):
    # bool, NaN, infinity, str and None are not numbers for the viewer.
    with pytest.raises(InvalidInputError, match="width"):
        check_number(value, "width")


def test_check_text_accepts_text():
    assert check_text("viridis", "name") == "viridis"


@pytest.mark.parametrize("value", ["", None, 3, b"a"])
def test_check_text_refuses(value):
    # The text must be a str and must not be empty.
    with pytest.raises(InvalidInputError, match="name"):
        check_text(value, "name")


def test_check_rgb_accepts_triple():
    assert check_rgb((0, 128, 255), "plain") == (0, 128, 255)
    assert check_rgb([1, 2, 3], "plain") == (1, 2, 3)


@pytest.mark.parametrize(
    "value",
    [(0, 0, 256), (-1, 0, 0), (1, 2), (1, 2, 3, 4), (0.5, 1, 2), (True, 0, 0), None, "abc"],
)
def test_check_rgb_refuses(value):
    # Three ints from 0 to 255 are the only valid value.
    with pytest.raises(InvalidInputError, match="plain"):
        check_rgb(value, "plain")


def test_check_codes_returns_list_of_int():
    result = check_codes([0, 1, 2**32 - 1], "codes")
    assert result == [0, 1, 2**32 - 1]
    assert isinstance(result, list)
    assert all(type(item) is int for item in result)


def test_check_codes_accepts_array():
    result = check_codes(np.array([3, 4], dtype=np.uint32), "codes")
    assert result == [3, 4]
    assert all(type(item) is int for item in result)


@pytest.mark.parametrize("values", [[True], [1.5], [-1], [2**32], [0, "1"], [None]])
def test_check_codes_refuses(values):
    # Each code must be an int from 0 to 2**32-1.
    with pytest.raises(InvalidInputError, match="codes"):
        check_codes(values, "codes")


def test_colour_settings_default_is_empty():
    settings = ColourSettings()
    assert settings.array is None
    assert settings.colour_map is None
    assert settings.range is None
    assert settings.scale is None
    assert settings.plain is None


def test_colour_settings_accepts_valid_values():
    settings = ColourSettings(
        array="a", colour_map="viridis", range=(0.0, 1.0), scale="log", plain=(1, 2, 3)
    )
    assert settings.array == "a"
    assert settings.scale == "log"


def test_colour_settings_is_frozen():
    settings = ColourSettings(array="a")
    with pytest.raises(dataclasses.FrozenInstanceError):
        settings.array = "b"


@pytest.mark.parametrize(
    ("kwargs", "field"),
    [
        ({"array": ""}, "array"),
        ({"array": 3}, "array"),
        ({"colour_map": ""}, "colour_map"),
        ({"range": (1.0, 1.0)}, "range"),
        ({"range": (2.0, 1.0)}, "range"),
        ({"range": (0.0,)}, "range"),
        ({"range": (0.0, float("nan"))}, "range"),
        ({"range": ("a", "b")}, "range"),
        ({"scale": "sqrt"}, "scale"),
        ({"plain": (0, 0, 256)}, "plain"),
        ({"plain": (0, 0)}, "plain"),
    ],
)
def test_colour_settings_refuses(kwargs, field):
    # Each error names the field that is not valid.
    with pytest.raises(InvalidInputError, match=field):
        ColourSettings(**kwargs)


def test_result_done_returns_none():
    assert raise_for_command_result({"status": "done"}, command="filter.edit") is None


def test_result_refused_raises_with_text():
    with pytest.raises(InvalidInputError) as info:
        raise_for_command_result({"status": "refused", "text": "x"}, command="filter.edit")
    assert "The viewer refused the command 'filter.edit': x" in str(info.value)


def test_result_unknown_status_is_protocol_error():
    with pytest.raises(ProtocolError) as info:
        raise_for_command_result({"status": "weird"}, command="filter.edit")
    assert "The viewer sent an unknown result status 'weird'." in str(info.value)


def test_parse_colour_map_names_gives_names():
    frame = _commands_frame([_map_entry(["viridis", "grey", "coolwarm"])])
    assert parse_colour_map_names(frame) == ("viridis", "grey", "coolwarm")


def test_parse_colour_map_names_without_enum_gives_none():
    entry = {
        "name": "colour.set_map",
        "input_schema": {"properties": {"name": {"type": "string"}}},
    }
    assert parse_colour_map_names(_commands_frame([entry])) is None


def test_parse_colour_map_names_without_entry_gives_none():
    other = {"name": "object.set_visible", "input_schema": {"properties": {}}}
    assert parse_colour_map_names(_commands_frame([other])) is None


@pytest.mark.parametrize("enum", ["viridis", [1, 2], ["a", 2], None])
def test_parse_colour_map_names_bad_enum_is_protocol_error(enum):
    # An enum that is not a list of str is a broken reply.
    with pytest.raises(ProtocolError):
        parse_colour_map_names(_commands_frame([_map_entry(enum)]))


def test_run_command_sends_request_and_returns_result():
    tool = StubTool(_result_frame({"status": "done"}))
    result = run_command(tool, "filter.edit", {"a": 1}, object_id="g1", timeout=5.0)
    assert result == {"status": "done"}
    ((build, kwargs),) = tool.calls
    assert kwargs["until"] == ("command.result",)
    assert kwargs["feature"] == FEATURE
    assert kwargs["object_id"] == "g1"
    assert kwargs["timeout"] == 5.0
    header = decode_frame(build(7)).header
    assert header["request"] == 7
    assert header["command"] == "filter.edit"
    assert header["input"] == {"a": 1}


def test_run_command_passes_cancel():
    tool = StubTool(_result_frame({"status": "done"}))
    cancel = object()
    run_command(tool, "filter.edit", {}, cancel=cancel)
    assert tool.calls[0][1]["cancel"] is cancel


def test_run_command_refused_raises():
    tool = StubTool(_result_frame({"status": "refused", "text": "no"}))
    with pytest.raises(InvalidInputError, match="no"):
        run_command(tool, "filter.edit", {})


@pytest.mark.parametrize("result", ["done", None, {"status": 3}, {"text": "x"}])
def test_run_command_bad_result_is_protocol_error(result):
    # The result must be a mapping with a str status.
    tool = StubTool(_result_frame(result))
    with pytest.raises(ProtocolError):
        run_command(tool, "filter.edit", {})


def test_list_colour_map_names_gives_names():
    tool = StubTool(_commands_frame([_map_entry(["viridis", "grey"])]))
    assert list_colour_map_names(tool) == ("viridis", "grey")
    ((build, kwargs),) = tool.calls
    assert kwargs["until"] == ("commands",)
    assert kwargs["feature"] == FEATURE
    assert decode_frame(build(3)).header["type"] == "commands.list"


# -- FilterResult ------------------------------------------------------------------------


def test_filter_result_counts():
    result = FilterResult.from_result({"status": "done", "passing": 100, "total": 1000})
    assert result == FilterResult(100, 1000, None)


def test_filter_result_empty():
    assert FilterResult.from_result({"status": "done"}) == FilterResult(None, None, None)


@pytest.mark.parametrize("state", ["unbuilt", "building"])
def test_filter_result_state(state):
    result = FilterResult.from_result({"status": "done", "state": state})
    assert result == FilterResult(None, None, state)


def test_filter_result_is_frozen():
    with pytest.raises(dataclasses.FrozenInstanceError):
        FilterResult().passing = 1  # type: ignore[misc]


@pytest.mark.parametrize("bad", [True, -1, 1.5, "10"])
def test_filter_result_bad_passing(bad):
    with pytest.raises(ProtocolError):
        FilterResult.from_result({"status": "done", "passing": bad})


@pytest.mark.parametrize("bad", [True, -1, 1.5, "10"])
def test_filter_result_bad_total(bad):
    with pytest.raises(ProtocolError):
        FilterResult.from_result({"status": "done", "total": bad})


@pytest.mark.parametrize("bad", [1, True, ["building"], b"building"])
def test_filter_result_bad_state(bad):
    with pytest.raises(ProtocolError):
        FilterResult.from_result({"status": "done", "state": bad})
