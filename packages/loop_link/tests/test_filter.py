"""Tests for the filter calls of ``Viewer`` (story 6)."""

from __future__ import annotations

import contextlib
import threading
import time
from collections.abc import Callable, Iterator
from typing import Any

import numpy as np
import pytest

from loop_link.commands import FilterResult
from loop_link.errors import BusyError, CancelledError, InvalidInputError, MissingFeatureError
from loop_link.objects import BlockModel, Categorical
from loop_link.testing import FakeViewer
from loop_link.viewer import Viewer

WAIT = 3.0
CELLS = (4, 4, 4)


def _grid() -> BlockModel:
    """Return a 4x4x4 block model with a number, a category and a text property."""
    grade = np.linspace(0.0, 1.0, 64)
    unit = Categorical(np.arange(64) % 3, ["a", "b", "c"])
    return BlockModel(
        (0, 0, 0),
        (1, 1, 1),
        CELLS,
        cell_properties={"grade": grade, "unit": unit},
        properties={"t": np.array(["a"] * 125)},
    )


def _buffers(fake: FakeViewer) -> dict[str, np.ndarray]:
    """Return a copy of the buffers of the object ``grid`` in the fake viewer."""
    return {name: array.copy() for name, array in fake.objects["grid"].buffers.items()}


def _runs(fake: FakeViewer) -> list[Any]:
    """Return the ``command.run`` frames that the fake viewer received."""
    return [f for f in fake.tool_requests if f.header.get("type") == "command.run"]


def _edits(fake: FakeViewer) -> list[dict[str, Any]]:
    """Return the edits of the ``filter.edit`` commands that the fake viewer received."""
    return [
        f.header["input"]["edit"] for f in _runs(fake) if f.header.get("command") == "filter.edit"
    ]


@pytest.fixture
def viewer(viewer_factory: Callable[..., Any]) -> Any:
    """Return a viewer that has sent the block model ``grid``."""
    made = viewer_factory(timeout=WAIT)
    made.send(_grid(), id="grid")
    made.flush()
    return made


# -- scenario 1: a range filter ----------------------------------------------------------


def test_filter_range_sets_filter(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_range("grid", "grade", 0.2, 0.8, keep_nan=True)
    state = fake_viewer.viewer.filters["grid"]["grade"]
    assert state["kind"] == "range"
    assert state["on"] is True
    assert state["minimum"] == 0.2
    assert state["maximum"] == 0.8
    assert state["keep_nan"] is True


def test_filter_range_edits_are_add_set_range_set_on(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_range("grid", "grade", 0.2, 0.8, keep_nan=True)
    edits = _edits(fake_viewer.viewer)
    assert [edit["kind"] for edit in edits] == ["add", "set_range", "set_on"]
    assert edits[1]["min"] == "0.2"
    assert edits[1]["max"] == "0.8"
    assert edits[1]["keep_nan"] is True
    assert edits[2]["on"] is True


# -- scenario 2: a values filter ---------------------------------------------------------


def test_filter_values_sets_filter(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_values("grid", "unit", [0, 2], exclude=True)
    state = fake_viewer.viewer.filters["grid"]["unit"]
    assert state["kind"] == "values"
    assert state["values"] == [0, 2]
    assert state["exclude"] is True
    assert state["on"] is True


def test_filter_values_edits_are_add_set_values_set_exclude_set_on(
    fake_viewer: Any, viewer: Any
) -> None:
    viewer.filter_values("grid", "unit", [0, 2], exclude=True)
    kinds = [edit["kind"] for edit in _edits(fake_viewer.viewer)]
    assert kinds == ["add", "set_values", "set_exclude", "set_on"]


# -- scenario 3: enable, remove, data does not change ------------------------------------


def test_filter_enable_false_turns_filter_off(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_range("grid", "grade", 0.2, 0.8)
    viewer.filter_enable("grid", "grade", False)
    assert fake_viewer.viewer.filters["grid"]["grade"]["on"] is False
    viewer.filter_enable("grid", "grade")
    assert fake_viewer.viewer.filters["grid"]["grade"]["on"] is True


def test_filter_remove_deletes_filter(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_range("grid", "grade", 0.2, 0.8)
    viewer.filter_remove("grid", "grade")
    assert "grade" not in fake_viewer.viewer.filters.get("grid", {})


def test_filter_calls_do_not_change_the_data(fake_viewer: Any, viewer: Any) -> None:
    fake = fake_viewer.viewer
    before = _buffers(fake)
    viewer.filter_range("grid", "grade", 0.2, 0.8)
    viewer.filter_enable("grid", "grade", False)
    viewer.filter_remove("grid", "grade")
    after = _buffers(fake)
    assert before.keys() == after.keys()
    for name, array in before.items():
        assert np.array_equal(array, after[name])


# -- scenario 4: bad input ---------------------------------------------------------------


@pytest.mark.parametrize("array", ["missing", "t"])
def test_filter_range_bad_array_names_the_array(viewer: Any, array: str) -> None:
    with pytest.raises(InvalidInputError) as info:
        viewer.filter_range("grid", array, 0.0, 1.0)
    assert array in str(info.value)


@pytest.mark.parametrize("array", ["missing", "t"])
def test_filter_values_bad_array_names_the_array(viewer: Any, array: str) -> None:
    with pytest.raises(InvalidInputError) as info:
        viewer.filter_values("grid", array, [0])
    assert array in str(info.value)


def test_filter_range_on_categorical_array_is_invalid(viewer: Any) -> None:
    with pytest.raises(InvalidInputError):
        viewer.filter_range("grid", "unit", 0.0, 1.0)


def test_filter_values_on_number_array_is_invalid(viewer: Any) -> None:
    with pytest.raises(InvalidInputError):
        viewer.filter_values("grid", "grade", [0])


def test_filter_range_minimum_above_maximum_sends_no_command(fake_viewer: Any, viewer: Any) -> None:
    with pytest.raises(InvalidInputError):
        viewer.filter_range("grid", "grade", 0.8, 0.2)
    assert _runs(fake_viewer.viewer) == []


@pytest.mark.parametrize("values", [[1.5], [True], [-1], [2**32]])
def test_filter_values_bad_value_is_invalid(viewer: Any, values: list[Any]) -> None:
    with pytest.raises(InvalidInputError):
        viewer.filter_values("grid", "unit", values)


# -- scenario 6: a new revision ----------------------------------------------------------


def test_filters_stay_after_a_new_revision(fake_viewer: Any, viewer: Any) -> None:
    viewer.filter_range("grid", "grade", 0.2, 0.8)
    before = fake_viewer.viewer.filters
    viewer.send(_grid(), id="grid")
    viewer.flush()
    assert fake_viewer.viewer.objects["grid"].header["revision"] == 2
    assert fake_viewer.viewer.filters == before


# -- refused add, cancel, viewer kinds, old protocol -------------------------------------


def test_filter_range_succeeds_when_add_is_refused(
    fake_viewer: Any, viewer: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = FakeViewer._filter_edit

    def refuse_add(self: FakeViewer, object_id: str, edit: dict[str, Any]) -> dict[str, Any]:
        if edit["kind"] == "add":
            return {"status": "refused", "text": "exists"}
        return original(self, object_id, edit)

    monkeypatch.setattr(FakeViewer, "_filter_edit", refuse_add)
    # The filter exists already, so the fake viewer needs it before the range is set.
    original(fake_viewer.viewer, "grid", {"kind": "add", "array": "grade"})
    viewer.filter_range("grid", "grade", 0.2, 0.8)
    state = fake_viewer.viewer.filters["grid"]["grade"]
    assert state["minimum"] == 0.2
    assert state["on"] is True


def test_cancel_filter_sends_one_cancel_edit(fake_viewer: Any, viewer: Any) -> None:
    viewer.cancel_filter("grid")
    runs = _runs(fake_viewer.viewer)
    assert len(runs) == 1
    assert runs[0].header["command"] == "filter.edit"
    assert runs[0].header["input"]["edit"] == {"kind": "cancel"}


def test_filter_works_with_headless_viewer() -> None:
    from loop_link.viewer import Viewer

    with FakeViewer(kind="headless") as fake:
        client = Viewer(port=fake.port, token="secret", timeout=WAIT)
        try:
            client.send(_grid(), id="grid")
            client.flush()
            client.filter_range("grid", "grade", 0.2, 0.8)
            assert fake.filters["grid"]["grade"]["on"] is True
        finally:
            client.close()


def test_filter_old_protocol_gives_missing_feature() -> None:
    from loop_link.viewer import Viewer

    with FakeViewer(protocol_minor=10) as fake:
        client = Viewer(port=fake.port, token="secret", timeout=WAIT)
        try:
            with pytest.raises(MissingFeatureError):
                client.filter_range("grid", "grade", 0.2, 0.8)
        finally:
            client.close()


# -- FR-014: the tool part works while another program holds the data connection ---------


def test_second_program_filters_without_data_connection(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    first = viewer_factory(timeout=WAIT)
    first.send(_grid(), id="grid")
    first.flush()
    second = viewer_factory(timeout=WAIT)
    with pytest.raises(BusyError):
        second.send(_grid(), id="other")
    second.filter_range("grid", "grade", 0.2, 0.8)
    assert fake_viewer.viewer.filters["grid"]["grade"]["on"] is True
    second.filter_remove("grid", "grade")
    assert "grade" not in fake_viewer.viewer.filters.get("grid", {})
    assert not second.connected


# -- the wait of the filter calls (protocol 1.20) -----------------------------------------

BIG = 1000


@pytest.fixture
def slow_viewer(
    tmp_path: Any, monkeypatch: pytest.MonkeyPatch
) -> Iterator[Callable[..., tuple[FakeViewer, Any]]]:
    """Return a function that starts a fake viewer and a viewer that has sent 1,000 cells."""
    opened: list[Any] = []

    def make(**options: Any) -> tuple[FakeViewer, Any]:
        fake = FakeViewer(**options)
        settings = tmp_path / f"geoviewer{len(opened)}"
        fake.write_instance_file(settings)
        monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
        viewer = Viewer(settings_dir=settings, heartbeat_interval=0.1, timeout=WAIT)
        opened.append((fake, viewer))
        grade = np.arange(BIG, dtype=np.float32)
        unit = Categorical(np.arange(BIG) % 3, ["a", "b", "c"])
        model = BlockModel(
            (0, 0, 0),
            (1, 1, 1),
            (10, 10, 10),
            cell_properties={"grade": grade, "unit": unit},
        )
        viewer.send(model, id="grid")
        viewer.flush()
        return fake, viewer

    yield make
    for fake, viewer in opened:
        with contextlib.suppress(Exception):  # errors on close do not matter at teardown
            viewer.close()
        fake.close()


def test_filter_wait_range_returns_counts_after_build(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4)
    started = time.monotonic()
    result = viewer.filter_range("grid", "grade", 0.0, 99.0)
    assert time.monotonic() - started >= 0.3
    assert result == FilterResult(100, 1000, None)


def test_filter_wait_values_returns_counts(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4)
    result = viewer.filter_values("grid", "unit", [0])
    assert result.total == 1000
    assert result.passing == 334
    assert result.state is None


def test_filter_wait_enable_and_remove_return_result(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.1)
    viewer.filter_range("grid", "grade", 0.0, 99.0)
    assert isinstance(viewer.filter_enable("grid", "grade", False), FilterResult)
    assert isinstance(viewer.filter_enable("grid", "grade", True), FilterResult)
    assert isinstance(viewer.filter_remove("grid", "grade"), FilterResult)


def test_filter_wait_progress_gets_fractions(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4)
    fractions: list[float] = []
    viewer.filter_range("grid", "grade", 0.0, 99.0, progress=fractions.append)
    assert fractions
    assert all(isinstance(value, float) for value in fractions)
    assert fractions[0] == 0.0
    assert fractions[-1] == 1.0


def test_filter_wait_false_returns_at_once(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4)
    started = time.monotonic()
    result = viewer.filter_range("grid", "grade", 0.0, 99.0, wait=False)
    assert time.monotonic() - started < 0.2
    assert result == FilterResult()
    # A next call works while the build goes on.
    assert isinstance(viewer.filter_enable("grid", "grade", True), FilterResult)


def test_filter_wait_cancel_event_raises_cancelled(slow_viewer: Any) -> None:
    fake, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4)
    cancel = threading.Event()
    timer = threading.Timer(0.1, cancel.set)
    timer.start()
    try:
        with pytest.raises(CancelledError):
            viewer.filter_range("grid", "grade", 0.0, 99.0, cancel=cancel)
    finally:
        timer.cancel()
    state = fake.filters.get("grid", {}).get("grade")
    assert state is None or state["on"] is False


def test_filter_wait_stop_gives_building_state(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=20, filter_delay=0.4, filter_stop_after=0.1)
    result = viewer.filter_range("grid", "grade", 0.0, 99.0)
    assert result.state == "building"


def test_filter_wait_old_viewer_gives_empty_result(slow_viewer: Any) -> None:
    _, viewer = slow_viewer(protocol_minor=18)
    assert viewer.filter_range("grid", "grade", 0.0, 99.0) == FilterResult()
