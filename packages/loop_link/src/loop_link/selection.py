"""The selection in the viewer: the ``Selection`` value and the ``SelectionWatcher``."""

from __future__ import annotations

import threading
import weakref
from collections.abc import Callable
from dataclasses import dataclass
from types import TracebackType
from typing import TYPE_CHECKING

import numpy as np

from loop_link.protocol import Frame, ProtocolError

if TYPE_CHECKING:
    from loop_link.viewer import Viewer

_KIND_FIELDS = {"cells": "cells", "nodes": "vertices"}


def _empty() -> np.ndarray:
    """Return an empty ``uint32`` array."""
    return np.empty(0, dtype=np.uint32)


@dataclass(frozen=True, eq=False)
class Selection:
    """The cells and vertices that are selected in the viewer.

    Parameters
    ----------
    object_id : str or None
        The identifier of the selected object. ``None`` when nothing is selected.
    revision : int or None
        The revision of the object when the selection was read.
    cells : numpy.ndarray
        The selected cell indices. A 1-D ``uint32`` array.
    vertices : numpy.ndarray
        The selected vertex indices. A 1-D ``uint32`` array.
    """

    object_id: str | None
    revision: int | None
    cells: np.ndarray
    vertices: np.ndarray

    @property
    def empty(self) -> bool:
        """Return True when no cell and no vertex is selected."""
        return self.cells.size == 0 and self.vertices.size == 0

    def same_as(self, other: Selection) -> bool:
        """Return True when ``other`` selects the same indices of the same object.

        The order of the indices and the revision are not compared.

        Parameters
        ----------
        other : Selection
            The selection to compare with.

        Returns
        -------
        bool
            True when the object ids are equal and the sorted cells and the sorted
            vertices are equal.
        """
        return (
            self.object_id == other.object_id
            and np.array_equal(np.sort(self.cells), np.sort(other.cells))
            and np.array_equal(np.sort(self.vertices), np.sort(other.vertices))
        )


def parse_selection_result(frame: Frame) -> Selection:
    """Return the selection in a ``selection.result`` frame.

    Parameters
    ----------
    frame : Frame
        The decoded ``selection.result`` frame.

    Returns
    -------
    Selection
        The selection. The entry kind ``cells`` gives ``cells``. The entry kind
        ``nodes`` gives ``vertices``.

    Raises
    ------
    ProtocolError
        If an entry has a count but no buffer, a buffer is not ``uint32``, or the
        length of a buffer is not the count.
    """
    header = frame.header
    arrays = {"cells": _empty(), "vertices": _empty()}
    for entry in header.get("entries") or []:
        field = _KIND_FIELDS.get(entry.get("kind"))
        if field is None:
            continue
        count = int(entry.get("count", 0))
        name = entry.get("indices")
        if name is None or name not in frame.buffers:
            if count:
                raise ProtocolError(
                    f"The selection entry '{entry.get('kind')}' has a count but no indices."
                )
            continue
        array = frame.buffers[name]
        if array.dtype != np.uint32:
            raise ProtocolError(
                f"The selection indices '{name}' are {array.dtype}. They must be uint32."
            )
        if array.ndim != 1 or array.shape[0] != count:
            raise ProtocolError(
                f"The selection indices '{name}' have shape {array.shape}. "
                f"The entry count is {count}."
            )
        arrays[field] = array
    return Selection(
        header.get("object"), header.get("revision"), arrays["cells"], arrays["vertices"]
    )


class SelectionWatcher:
    """Poll the selection of a viewer and call a function when it changes.

    Only ``Viewer.watch_selection`` makes a watcher. The thread holds only a weak
    reference to the viewer.

    Parameters
    ----------
    viewer_ref : weakref.ref
        A weak reference to the viewer.
    callback : callable
        The function to call with the new ``Selection`` at each change.
    interval : float
        The seconds between two polls.
    last : Selection
        The selection that counts as already seen.
    """

    def __init__(
        self,
        viewer_ref: weakref.ref[Viewer],
        callback: Callable[[Selection], None],
        interval: float,
        last: Selection,
    ) -> None:
        self._viewer_ref = viewer_ref
        self._callback = callback
        self._interval = interval
        self._last = last
        self._stop_event = threading.Event()
        self._error: BaseException | None = None
        self._thread = threading.Thread(
            target=self._run, name="loop-link-selection-poll", daemon=True
        )

    def start(self) -> None:
        """Start the poll thread."""
        self._thread.start()

    @property
    def running(self) -> bool:
        """Return True when the poll thread is alive."""
        return self._thread.is_alive()

    @property
    def error(self) -> BaseException | None:
        """Return the exception that ended the thread, or None."""
        return self._error

    def stop(self) -> None:
        """Stop the poll thread and wait for it. A second call does no harm."""
        self._stop_event.set()
        if self._thread.is_alive() and self._thread is not threading.current_thread():
            self._thread.join(2 * self._interval + 1)

    def __enter__(self) -> SelectionWatcher:  # noqa: PYI034
        """Return the watcher."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Stop the watcher."""
        self.stop()

    def _run(self) -> None:
        """Poll until stopped, until the viewer is gone, or until an error."""
        while not self._stop_event.wait(self._interval):
            viewer = self._viewer_ref()
            if viewer is None:
                return
            try:
                selection = viewer.selection()
            except Exception as error:  # noqa: BLE001
                self._error = error
                return
            finally:
                del viewer
            if selection.same_as(self._last):
                continue
            self._last = selection
            try:
                self._callback(selection)
            except Exception as error:  # noqa: BLE001
                self._error = error
                return
