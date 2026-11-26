"""Tests for the memory accounting of the fault evaluation (D4, D5, FR-010, T022, T047).

* The workspace of the digest of the points is bounded. It is far smaller than the
  arrays that are hashed.
* No restoration stays after a call.
* One call keeps each batch of the caller once and at most one frame-step batch for
  each fault set and size (FR-008, FR-010, C-018).

The large case with 10 million points and its tracemalloc peak is not part of normal pytest
(T029, T030).
"""

from __future__ import annotations

import gc
import tracemalloc
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from LoopStructural.modelling.features import _evaluation
from tests.fault_evaluation.models import build_model, sample_points

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


# ---------------------------------------------------------------------------
# Held restorations during a call (FR-010)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["evaluate_model", "evaluate_model_gradient"])
def test_scope_holds_at_most_one_restoration_per_fault_set(monkeypatch, method):
    """The scope holds at most one restored batch for each fault set and batch size.

    The fault sets include the fault sets of the fault frames, which the frame
    steps restore through older faults. Per OP1 the memory is counted per actual
    restored point, so a fault set may hold one N-point batch and one batch of
    4N tetrahedron corners. Before the fix the scope kept every frame step batch.
    """
    model = build_model("cutting_faults")
    points = sample_points("cutting_faults", 5000, seed=1)
    largest: dict = {}
    held = []
    restore = _evaluation.EvaluationContext.restore

    def spy(self, feature, batch, reverse=False):
        restored = restore(self, feature, batch, reverse)
        if feature.faults_enabled and feature.faults:
            batch = np.asarray(batch)
            slot = (tuple(map(id, feature.faults)), bool(reverse), batch.shape)
            largest[slot] = max(largest.get(slot, 0), batch.nbytes)
        held.append(sum(array.nbytes for array in self._restorations.values()))
        return restored

    monkeypatch.setattr(_evaluation.EvaluationContext, "restore", spy)
    getattr(model, method)(points.copy())
    # cutting_faults: the model fault set and the frame fault set of the younger
    # fault. The frame set has N points and 4N corners for each frame step.
    fault_sets = {slot[:2] for slot in largest}
    assert len(fault_sets) == 2
    assert len(largest) == 3
    assert max(held) <= sum(largest.values())


def test_caller_batches_of_one_size_do_not_replace_each_other():
    """Two caller batches of the same size on one fault set are each restored once (FR-008).

    Review round 2: with one slot for each fault set and size, the calls e, e + 1, e
    applied the fault three times.
    """
    from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation

    model = build_model("cutting_faults")
    member = model.get_feature_by_name("fault_1")[0]
    points = sample_points("cutting_faults", 500, seed=1)
    with EvaluationInstrumentation() as counter, _evaluation.evaluation_scope(model):
        member.evaluate_value(points)
        member.evaluate_value(points + 1.0)
        member.evaluate_value(points)
    assert counter.snapshot()["fault_calls_total"] == 2 * len(member.faults)
