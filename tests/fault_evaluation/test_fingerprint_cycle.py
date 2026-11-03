"""Regression test: a cycle back to the root must not hide the deep state of the root.

The root of ``dependency_fingerprint`` is first written shallow. When a region of a fault
reads the root, the root is met again as a dependency. That visit must write the full
state, so a change of the field of the root must change the digest.
"""

from __future__ import annotations

import threading

from LoopStructural.modelling.features._evaluation_cache import dependency_fingerprint
from LoopStructural.modelling.features._region import Region
from tests.fault_evaluation.models import build_model


def _fingerprint_with_timeout(feature) -> bytes | None:
    """Call the function in a thread. Fail when it does not end."""
    result: list = []
    thread = threading.Thread(target=lambda: result.append(dependency_fingerprint(feature)))
    thread.daemon = True
    thread.start()
    thread.join(60)
    assert not thread.is_alive(), "dependency_fingerprint did not end (cycle?)"
    return result[0]


def test_root_state_is_hashed_when_a_region_reads_the_root():
    model = build_model("f1_u0")
    t = model.get_feature_by_name("group_0")
    t.faults[0].regions.append(Region(t, 0.0, 1))

    first = _fingerprint_with_timeout(t)
    assert first is not None
    # No change: the digest is the same and the call ends.
    assert _fingerprint_with_timeout(t) == first

    t.interpolator.c[:] += 1
    assert _fingerprint_with_timeout(t) != first
