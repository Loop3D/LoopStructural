"""Speed checks for the loop_link send path."""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

from loop_link.objects import BlockModel
from loop_link.testing import FakeViewer

POLL_INTERVAL = 0.01  # seconds between inventory checks
GIVE_UP_AFTER = 10.0  # seconds before the wait is a failure
SPEED_LIMIT = 5.0  # seconds, SC-001


def _stored_grade(viewer: FakeViewer) -> np.ndarray:
    """Return the grade buffer of the object with id 'grid'."""
    stored = viewer.objects["grid"]
    for entry in stored.header["data"]:
        if entry.get("name") == "grade":
            return stored.buffers[entry["buffer"]]
    raise AssertionError("the stored object has no grade data")


@pytest.mark.slow
def test_send_million_cells_under_five_seconds(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """SC-001: a 1,000,000-cell block model reaches the viewer in under 5 s."""
    grade = np.random.default_rng(0).random(1_000_000, dtype=np.float32)
    block = BlockModel(
        origin=(0, 0, 0),
        cell_sizes=(1, 1, 1),
        cells=(100, 100, 100),
        cell_properties={"grade": grade},
    )
    viewer = viewer_factory()

    start = time.perf_counter()
    viewer.send(block, id="grid")
    deadline = start + GIVE_UP_AFTER
    while not any(item["id"] == "grid" for item in fake_viewer.viewer.inventory()):
        if time.perf_counter() > deadline:
            pytest.fail("the viewer did not list 'grid' within 10 s")
        time.sleep(POLL_INTERVAL)
    elapsed = time.perf_counter() - start
    print(f"\nSC-001 send time: {elapsed:.3f} s")

    assert elapsed < SPEED_LIMIT
    assert np.array_equal(_stored_grade(fake_viewer.viewer), grade)
