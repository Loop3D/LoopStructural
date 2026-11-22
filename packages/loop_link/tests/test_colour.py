"""Tests for the colour calls of ``Viewer`` (story 7)."""

from __future__ import annotations

import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from loop_link.commands import ColourSettings
from loop_link.errors import InvalidInputError, MissingFeatureError, ViewerRefusedError
from loop_link.objects import BlockModel, Categorical
from loop_link.testing import FakeViewer

WAIT = 3.0
COLOUR_COMMANDS = (
    "colour.set_active",
    "colour.set_map",
    "colour.set_range",
    "colour.set_scale",
    "object.set_colour",
)


def _grid() -> BlockModel:
    """Return a 4x4x4 block model with a number, a category and a boolean cell property."""
    count = 4 * 4 * 4
    return BlockModel(
        (0, 0, 0),
        (1, 1, 1),
        (4, 4, 4),
        cell_properties={
            "grade": np.linspace(0.0, 1.0, count),
            "unit": Categorical(np.arange(count) % 2, ["a", "b"]),
            "flag": np.arange(count) % 2 == 0,
        },
    )


def _commands(fake: FakeViewer, *names: str) -> list[dict[str, Any]]:
    """Return the ``command.run`` frames that the fake received, as header dicts.

    If ``names`` is given, return only the frames of these commands.
    """
    headers = [
        frame.header for frame in fake.tool_requests if frame.header.get("type") == "command.run"
    ]
    if names:
        headers = [header for header in headers if header.get("command") in names]
    return headers


def _make(
    tmp_path: Path, viewer_factory: Callable[..., Any], **options: Any
) -> tuple[FakeViewer, Any]:
    """Return a second fake viewer and a viewer that is connected to it.

    The viewer has sent the block model ``grid``.
    """
    fake = FakeViewer(**options)
    fake.write_instance_file(tmp_path / "other")
    viewer = viewer_factory(settings_dir=tmp_path / "other", timeout=WAIT)
    viewer.send(_grid(), id="grid")
    viewer.flush()
    return fake, viewer


@pytest.fixture
def viewer(viewer_factory: Callable[..., Any]) -> Any:
    """Return a viewer that has sent the block model ``grid``."""
    made = viewer_factory(timeout=WAIT)
    made.send(_grid(), id="grid")
    made.flush()
    return made


def test_set_colour_property_sets_array(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 1: the array name goes to the viewer."""
    viewer.set_colour_property("grid", "grade")
    assert fake_viewer.viewer.colouring["grid"]["array"] == "grade"
    sent = _commands(fake_viewer.viewer, "colour.set_active")
    assert sent[-1]["input"] == {"id": "grid", "array": "grade"}


def test_set_colour_property_none_leaves_array_out(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 1: ``array=None`` sends no ``array`` key and the viewer shows no array."""
    viewer.set_colour_property("grid", "grade")
    viewer.set_colour_property("grid", None)
    sent = _commands(fake_viewer.viewer, "colour.set_active")
    assert "array" not in sent[-1]["input"]
    assert fake_viewer.viewer.colouring["grid"]["array"] is None


def test_set_colour_map_and_names(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 2: the colour map name goes to the viewer. The names come from the viewer."""
    assert viewer.colour_map_names() == ("viridis", "grey", "coolwarm")
    viewer.set_colour_map("grid", "coolwarm")
    assert fake_viewer.viewer.colouring["grid"]["colour_map"] == "coolwarm"


def test_set_colour_range_sends_text(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 3: the limits go as text."""
    viewer.set_colour_property("grid", "grade")
    viewer.set_colour_range("grid", 0.1, 0.9)
    sent = _commands(fake_viewer.viewer, "colour.set_range")
    assert sent[-1]["input"] == {"id": "grid", "min": "0.1", "max": "0.9"}
    assert fake_viewer.viewer.colouring["grid"]["range"] == (0.1, 0.9)


def test_clear_colour_range(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 3: ``clear_colour_range`` removes the range."""
    viewer.set_colour_property("grid", "grade")
    viewer.set_colour_range("grid", 0.1, 0.9)
    viewer.clear_colour_range("grid")
    assert fake_viewer.viewer.colouring["grid"]["range"] is None
    assert _commands(fake_viewer.viewer, "colour.clear_range")


@pytest.mark.parametrize("limits", [(0.9, 0.1), (0.5, 0.5)])
def test_set_colour_range_refuses_bad_limits(
    viewer: Any, fake_viewer: Any, limits: tuple[float, float]
) -> None:
    """Scenario 3: a low limit that is not less than the high limit sends no frame."""
    with pytest.raises(InvalidInputError):
        viewer.set_colour_range("grid", *limits)
    assert not _commands(fake_viewer.viewer, "colour.set_range")


def test_set_colour_scale(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 4: the scale goes to the viewer. The input key is ``object``."""
    viewer.set_colour_property("grid", "grade")
    viewer.set_colour_scale("grid", "log")
    sent = _commands(fake_viewer.viewer, "colour.set_scale")
    assert sent[-1]["input"] == {"object": "grid", "scale": "log"}
    assert fake_viewer.viewer.colouring["grid"]["scale"] == "log"


def test_set_colour_scale_refuses_unknown_scale(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 4: an unknown scale sends no frame."""
    with pytest.raises(InvalidInputError):
        viewer.set_colour_scale("grid", "sqrt")
    assert not _commands(fake_viewer.viewer, "colour.set_scale")


def test_set_object_colour(viewer: Any, fake_viewer: Any) -> None:
    """The plain colour goes as three integers. ``None`` goes as null."""
    viewer.set_object_colour("grid", (255, 0, 0))
    sent = _commands(fake_viewer.viewer, "object.set_colour")
    assert sent[-1]["input"] == {"id": "grid", "colour": [255, 0, 0]}
    assert fake_viewer.viewer.colouring["grid"]["plain"] == [255, 0, 0]
    viewer.set_object_colour("grid", None)
    sent = _commands(fake_viewer.viewer, "object.set_colour")
    assert sent[-1]["input"] == {"id": "grid", "colour": None}
    assert fake_viewer.viewer.colouring["grid"]["plain"] is None


def test_set_object_colour_refuses_bad_value(viewer: Any, fake_viewer: Any) -> None:
    """A value above 255 sends no frame."""
    with pytest.raises(InvalidInputError):
        viewer.set_object_colour("grid", (256, 0, 0))
    assert not _commands(fake_viewer.viewer, "object.set_colour")


@pytest.mark.parametrize("array", ["unit", "flag"])
def test_fixed_colouring_refuses_map_range_scale(viewer: Any, array: str) -> None:
    """Scenario 5: an array with a fixed colouring refuses a map, a range and a scale."""
    viewer.set_colour_property("grid", array)
    with pytest.raises(InvalidInputError, match="fixed colouring"):
        viewer.set_colour_map("grid", "viridis")
    with pytest.raises(InvalidInputError, match="fixed colouring"):
        viewer.set_colour_range("grid", 0.0, 1.0)
    with pytest.raises(InvalidInputError, match="fixed colouring"):
        viewer.set_colour_scale("grid", "log")


def test_unknown_colour_map_lists_names(viewer: Any, fake_viewer: Any) -> None:
    """Scenario 6: an unknown name gives an error with the known names and sends no frame."""
    with pytest.raises(InvalidInputError) as caught:
        viewer.set_colour_map("grid", "nope")
    for name in ("viridis", "grey", "coolwarm"):
        assert name in str(caught.value)
    assert not _commands(fake_viewer.viewer, "colour.set_map")


def test_old_viewer_decides_colour_map_name(
    fake_viewer: Any, tmp_path: Path, viewer_factory: Callable[..., Any]
) -> None:
    """Scenario 6: a viewer of minor version 17 has no names. It gets the name and decides."""
    fake, viewer = _make(tmp_path, viewer_factory, protocol_minor=17)
    try:
        viewer.set_colour_map("grid", "nope")
        assert fake.colouring["grid"]["colour_map"] == "nope"
        with pytest.raises(MissingFeatureError):
            viewer.colour_map_names()
    finally:
        fake.close()


def test_apply_colour_sets_all_fields_in_order(viewer: Any, fake_viewer: Any) -> None:
    """``apply_colour`` sends the five commands in a fixed order."""
    settings = ColourSettings(
        array="grade", colour_map="viridis", range=(0.0, 1.0), scale="linear", plain=(1, 2, 3)
    )
    viewer.apply_colour("grid", settings)
    sent = _commands(fake_viewer.viewer, *COLOUR_COMMANDS)
    assert tuple(header["command"] for header in sent) == COLOUR_COMMANDS
    state = fake_viewer.viewer.colouring["grid"]
    assert state["array"] == "grade"
    assert state["colour_map"] == "viridis"
    assert state["range"] == (0.0, 1.0)
    assert state["scale"] == "linear"
    assert state["plain"] == [1, 2, 3]


def test_apply_colour_sends_only_given_fields(viewer: Any, fake_viewer: Any) -> None:
    """A ``ColourSettings`` with only ``scale`` sends one command."""
    viewer.apply_colour("grid", ColourSettings(scale="log"))
    sent = _commands(fake_viewer.viewer, *COLOUR_COMMANDS)
    assert [header["command"] for header in sent] == ["colour.set_scale"]


def test_colour_settings_refuses_bad_field_before_any_frame(viewer: Any, fake_viewer: Any) -> None:
    """A bad field of ``ColourSettings`` gives an error before any frame goes out."""
    with pytest.raises(InvalidInputError):
        viewer.apply_colour("grid", ColourSettings(scale="sqrt"))
    assert not _commands(fake_viewer.viewer, *COLOUR_COMMANDS)


def test_headless_viewer_has_no_colour_scale(
    fake_viewer: Any, tmp_path: Path, viewer_factory: Callable[..., Any]
) -> None:
    """A viewer without a window refuses the colour scale. The other colour calls work."""
    fake, viewer = _make(tmp_path, viewer_factory, kind="headless")
    try:
        with pytest.raises(MissingFeatureError):
            viewer.set_colour_scale("grid", "log")
        viewer.set_colour_property("grid", "grade")
        viewer.set_colour_map("grid", "grey")
        viewer.set_colour_range("grid", 0.0, 1.0)
        viewer.clear_colour_range("grid")
        viewer.set_object_colour("grid", (1, 2, 3))
        state = fake.colouring["grid"]
        assert state["array"] == "grade"
        assert state["colour_map"] == "grey"
        assert state["plain"] == [1, 2, 3]
    finally:
        fake.close()


def test_unknown_object_is_refused(viewer: Any) -> None:
    """A colour call for an object that is not in the scene gives ``InvalidInputError``."""
    with pytest.raises(InvalidInputError):
        viewer.set_colour_property("missing", "grade")


def _fresh_viewer(viewer_factory: Callable[..., Any], **options: Any) -> Any:
    """Return a viewer that has not sent anything. ``options`` go to the viewer."""
    return viewer_factory(timeout=WAIT, **options)


def test_send_colour_applies_all_settings(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """``send`` with ``colour`` sends the object and then applies all five settings."""
    viewer = _fresh_viewer(viewer_factory)
    settings = ColourSettings(
        array="grade", colour_map="coolwarm", range=(0.0, 1.0), scale="log", plain=(9, 8, 7)
    )
    revision = viewer.send(_grid(), id="grid", colour=settings)
    assert revision == 1
    assert "grid" in fake_viewer.viewer.objects
    state = fake_viewer.viewer.colouring["grid"]
    assert state["array"] == "grade"
    assert state["colour_map"] == "coolwarm"
    assert state["range"] == (0.0, 1.0)
    assert state["scale"] == "log"
    assert state["plain"] == [9, 8, 7]


def test_send_colour_commands_follow_apply_colour_order(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """The colour commands of ``send`` arrive in the order of ``apply_colour``."""
    viewer = _fresh_viewer(viewer_factory)
    settings = ColourSettings(
        array="grade", colour_map="viridis", range=(0.0, 1.0), scale="linear", plain=(1, 2, 3)
    )
    viewer.send(_grid(), id="grid", colour=settings)
    sent = _commands(fake_viewer.viewer, *COLOUR_COMMANDS)
    assert tuple(header["command"] for header in sent) == COLOUR_COMMANDS


def test_send_colour_object_is_in_viewer_before_first_colour_command(
    fake_viewer: Any, viewer_factory: Callable[..., Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The viewer holds the object when the first colour command arrives."""
    fake = fake_viewer.viewer
    seen: list[bool] = []
    original = FakeViewer._tool_reply

    def record(self: FakeViewer, frame: Any) -> list[Any]:
        if (
            frame.header.get("type") == "command.run"
            and frame.header.get("command") == "colour.set_active"
        ):
            seen.append("grid" in self.objects)
        return original(self, frame)

    monkeypatch.setattr(FakeViewer, "_tool_reply", record)
    viewer = _fresh_viewer(viewer_factory)
    viewer.send(_grid(), id="grid", colour=ColourSettings(array="grade"))
    assert seen == [True]
    assert fake.colouring["grid"]["array"] == "grade"


def test_send_colour_failed_command_keeps_object(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """A colour command that fails raises an error. The object stays in the viewer."""
    viewer = _fresh_viewer(viewer_factory)
    with pytest.raises(InvalidInputError):
        viewer.send(_grid(), id="grid", colour=ColourSettings(array="missing"))
    assert "grid" in fake_viewer.viewer.objects
    assert viewer.revisions["grid"] == 1


def test_send_colour_refused_object_sends_no_colour_command(
    tmp_path: Path, viewer_factory: Callable[..., Any]
) -> None:
    """A refused object raises ``ViewerRefusedError`` and no colour command goes out."""
    fake = FakeViewer(reply_errors={"grid": "bad"})
    fake.write_instance_file(tmp_path / "other")
    viewer = viewer_factory(settings_dir=tmp_path / "other", timeout=WAIT)
    try:
        with pytest.raises(ViewerRefusedError):
            viewer.send(_grid(), id="grid", colour=ColourSettings(array="grade"))
        assert not _commands(fake, *COLOUR_COMMANDS)
    finally:
        fake.close()


def test_send_colour_applies_when_another_object_was_refused(
    tmp_path: Path, viewer_factory: Callable[..., Any]
) -> None:
    """The viewer refuses ``other``. The colour of the accepted ``grid`` is still applied."""
    fake = FakeViewer(reply_errors={"other": "bad"})
    fake.write_instance_file(tmp_path / "other")
    viewer = viewer_factory(settings_dir=tmp_path / "other", timeout=WAIT)
    try:
        viewer.send(_grid(), id="other")  # no flush: the viewer refuses it later
        try:
            viewer.send(_grid(), id="grid", colour=ColourSettings(array="grade"))
        except ViewerRefusedError as error:
            # A refusal about "other" is allowed after the colour is applied.
            assert all(item.object_id == "other" for item in error.refusals)
        assert fake.colouring["grid"]["array"] == "grade"
    finally:
        fake.close()


def test_send_without_colour_sends_no_tool_frame(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """``colour=None`` sends no tool frame and opens no tool connection."""
    viewer = _fresh_viewer(viewer_factory)
    viewer.send(_grid(), id="grid", colour=None)
    viewer.flush()
    assert fake_viewer.viewer.tool_requests == []
    assert not viewer.tool_connected


def test_send_colour_wrong_type_is_refused_before_send(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """A ``colour`` that is not a ``ColourSettings`` gives an error. No object goes out."""
    viewer = _fresh_viewer(viewer_factory)
    with pytest.raises(InvalidInputError):
        viewer.send(_grid(), id="grid", colour={"array": "grade"})
    viewer.flush()
    assert not fake_viewer.viewer.objects


SPEED_LIMIT = 0.250  # seconds for one filter or colour call, SC-007


@pytest.mark.slow
def test_sc007_speed_million_cells(fake_viewer: Any, viewer_factory: Callable[..., Any]) -> None:
    """SC-007: each filter or colour call on 1 000 000 cells takes less than 0.250 s."""
    viewer = viewer_factory(timeout=30.0)
    model = BlockModel(
        origin=(0, 0, 0),
        cell_sizes=(1, 1, 1),
        cells=(100, 100, 100),
        cell_properties={"grade": np.random.default_rng(0).random(1_000_000, dtype=np.float32)},
    )
    viewer.send(model, id="grid")
    viewer.flush()
    viewer.set_colour_property("grid", "grade")  # warm-up call, not timed
    calls = {
        "set_colour_property": lambda: viewer.set_colour_property("grid", "grade"),
        "set_colour_map": lambda: viewer.set_colour_map("grid", "viridis"),
        "set_colour_range": lambda: viewer.set_colour_range("grid", 0.0, 1.0),
        "clear_colour_range": lambda: viewer.clear_colour_range("grid"),
        "set_colour_scale": lambda: viewer.set_colour_scale("grid", "linear"),
        "set_object_colour": lambda: viewer.set_object_colour("grid", (1, 2, 3)),
        "filter_range": lambda: viewer.filter_range("grid", "grade", 0.5, 1.0),
        "filter_enable": lambda: viewer.filter_enable("grid", "grade", False),
        "filter_remove": lambda: viewer.filter_remove("grid", "grade"),
    }
    for name, call in calls.items():
        start = time.perf_counter()
        call()
        elapsed = time.perf_counter() - start
        assert elapsed < SPEED_LIMIT, f"{name} took {elapsed:.3f} s"
