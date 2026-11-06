"""Push many objects to GeoViewer in one call.

The module has no Qt code and no module state. The progress callback runs on the
calling thread.
"""

from __future__ import annotations

import threading
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from loop_link.errors import InvalidInputError, ObjectTooLargeError, Refusal

if TYPE_CHECKING:
    from loop_link.viewer import Viewer


@dataclass
class PushResult:
    """Result of one call to `push_objects`.

    Attributes
    ----------
    sent : int
        Number of objects sent to the viewer.
    removed : int
        Number of old objects removed from the viewer.
    refused : list of Refusal
        Messages that the viewer refused.
    skipped : list of tuple of str
        Pairs of object id and reason for each object that was not sent.
    cancelled : bool
        True when the cancel event stopped the push.
    """

    sent: int = 0
    removed: int = 0
    refused: list[Refusal] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    cancelled: bool = False

    def summary(self) -> str:
        """Return a short text that describes the result.

        Returns
        -------
        str
            The summary text.
        """
        text = f"{self.sent} objects sent to GeoViewer"
        if self.removed:
            text += f", {self.removed} old objects removed"
        if self.skipped:
            text += f", {len(self.skipped)} objects skipped"
        text += "."
        if self.refused:
            text += f" The viewer refused {len(self.refused)} messages."
        return text


def push_objects(
    viewer: Viewer,
    objects: Iterable[tuple[str, Any]],
    *,
    replace_prefix: str | None = None,
    progress: Callable[[str], None] | None = None,
    cancel: threading.Event | None = None,
) -> PushResult:
    """Send objects to the viewer and report the result.

    Parameters
    ----------
    viewer : Viewer
        The viewer connection. The function connects it when it is not connected.
    objects : iterable of tuple of str and object
        Pairs of object id and object.
    replace_prefix : str, optional
        When given, the function removes each old object in the viewer inventory
        whose id starts with this prefix and that this call did not send.
    progress : callable, optional
        Called with a text before each object is sent. It runs on the calling thread.
    cancel : threading.Event, optional
        When this event is set, the push stops before the next object.

    Returns
    -------
    PushResult
        Counts, skipped objects and refusals.

    Raises
    ------
    ViewerClosedError
        If the viewer closes during the push.
    BusyError
        If the viewer is busy.
    """
    result = PushResult()
    if not viewer.connected:
        viewer.connect()
    sent_ids: set[str] = set()
    for object_id, obj in objects:
        if cancel is not None and cancel.is_set():
            result.cancelled = True
            break
        if progress is not None:
            progress(f"Sending {object_id}...")
        try:
            viewer.send(obj, id=object_id)
        except (InvalidInputError, ObjectTooLargeError) as error:
            result.skipped.append((object_id, str(error)))
            continue
        sent_ids.add(object_id)
        result.sent += 1
    if replace_prefix:
        for entry in viewer.inventory():
            if entry.id.startswith(replace_prefix) and entry.id not in sent_ids:
                viewer.remove(entry.id)
                result.removed += 1
    result.refused = viewer.flush(raise_on_refusal=False)
    return result
