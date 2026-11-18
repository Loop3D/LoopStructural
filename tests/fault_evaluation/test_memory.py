"""Tests for the memory accounting of the fault evaluation (D4, D5, FR-010, T022, T047).

* The workspace of the digest of the points is bounded. It is far smaller than the
  arrays that are hashed.
* No restoration stays after a call.

The large case with 10 million points and the peak RSS is not part of normal pytest
(T029, T030).
"""

from __future__ import annotations

import gc
import tracemalloc
from collections.abc import Callable
from pathlib import Path

import numpy as np

from LoopStructural.modelling.features import _evaluation
from tests.fault_evaluation.models import build_model

MEGABYTE = 1024 * 1024
DATA = Path(__file__).parent / "data"


def _peak(call: Callable[[], object]) -> int:
    """Return the peak of the Python and NumPy allocations of ``call`` in bytes."""
    call()  # warm call: imports and lazy tables are not counted
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        call()
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak


# ---------------------------------------------------------------------------
# Bounded hashing workspace
# ---------------------------------------------------------------------------


def test_points_digest_uses_a_bounded_workspace():
    """The digest of 24 MB of points allocates far less than the points.

    A copy of the array (for example ``tobytes``) would need 24 MB. The limit of 6 MB
    leaves room for the chunk buffer, the digest state and Python objects.
    """
    points = np.linspace(0.0, 1.0, 3_000_000).reshape(-1, 3)
    size = points.nbytes
    assert size >= 24_000_000
    peak = _peak(lambda: _evaluation._digest(points))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"


def test_points_digest_of_a_view_uses_a_bounded_workspace():
    """The digest of a large non-contiguous view must not copy the whole array.

    A Fortran-order array is a plausible input (for example a transposed array). Chunks
    of bounded size must be copied, not the full array.
    """
    points = np.asfortranarray(np.linspace(0.0, 1.0, 3_000_000).reshape(-1, 3))
    assert not points.flags.c_contiguous
    size = points.nbytes
    peak = _peak(lambda: _evaluation._digest(points))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"


def test_points_digest_is_the_same_for_the_chunked_view():
    """The chunked digest of a non-contiguous array equals the digest of its C-order copy."""
    contiguous = np.linspace(0.0, 1.0, 600_000).reshape(-1, 3)
    for variant in (np.asfortranarray(contiguous), np.repeat(contiguous, 2, axis=0)[::2]):
        assert _evaluation._digest(variant) == _evaluation._digest(contiguous)


# ---------------------------------------------------------------------------
# No restoration after a call
# ---------------------------------------------------------------------------


def test_no_restoration_stays_after_the_call():
    """After ``evaluate_model``, no scope is active and no context object is alive."""
    model = build_model("f1_u0")
    with np.load(DATA / "f1_u0.npz", allow_pickle=False) as archive:
        points = archive["points"].copy()
    model.evaluate_model(points)
    assert _evaluation._ACTIVE.get() is None
    gc.collect()
    contexts = [o for o in gc.get_objects() if isinstance(o, _evaluation.EvaluationContext)]
    assert contexts == [], "an EvaluationContext stays alive after the call"
