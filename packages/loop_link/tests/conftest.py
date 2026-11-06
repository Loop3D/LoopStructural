"""Shared test support for the loop_link tests."""

from __future__ import annotations

import contextlib
import os
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from loop_link.testing import FakeViewer


@dataclass
class FakeViewerHandle:
    """A running fake viewer and the settings folder that holds its instance file."""

    viewer: FakeViewer
    settings: Path


@pytest.fixture
def fake_viewer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[FakeViewerHandle]:
    """Start a fake viewer and publish its instance file in a temporary settings folder."""
    viewer = FakeViewer()
    settings = tmp_path / "geoviewer"
    try:
        viewer.write_instance_file(settings)
        monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
        yield FakeViewerHandle(viewer=viewer, settings=settings)
    finally:
        viewer.close()


@pytest.fixture
def viewer_factory(fake_viewer: FakeViewerHandle) -> Iterator[Callable[..., Any]]:
    """Return a function that makes viewers. Each made viewer is closed at teardown."""

    from loop_link.viewer import Viewer  # Viewer is built by a later task, so import it here

    made: list[Any] = []

    def make(**kwargs: Any) -> Any:
        options = {"settings_dir": fake_viewer.settings, "heartbeat_interval": 0.1}
        options.update(kwargs)
        viewer = Viewer(**options)
        made.append(viewer)
        return viewer

    yield make
    for viewer in made:
        with contextlib.suppress(Exception):  # errors on close do not matter at teardown
            viewer.close()


def geoviewer_fixtures() -> Path | None:
    """Return the folder in LOOP_LINK_GEOVIEWER_FIXTURES, or None when it is not a folder."""
    value = os.environ.get("LOOP_LINK_GEOVIEWER_FIXTURES")
    if not value:
        return None
    path = Path(value)
    return path if path.is_dir() else None
