"""Fixtures for the fault evaluation tests."""

from __future__ import annotations

import pytest

from LoopStructural import LoopStructuralConfig

PYTHON_PATH_MODULES = (
    "test_context",
    "test_gradient_context",
    "test_memory",
    "test_concurrency",
    "test_no_retention",
    "test_fault_frame_lookup",
)


@pytest.fixture(autouse=True)
def _python_path_for_internals(request, monkeypatch):
    """Turn off the native fault kernel for the internal path tests.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The request of the test that uses this fixture.
    monkeypatch : pytest.MonkeyPatch
        Patches the configuration for the duration of the test.

    Notes
    -----
    These tests count or patch calls of the Python path
    (``FaultSegment.apply_to_points`` or the grid lookups). The native
    path does not make these calls. See spec 007 plan, Answers, T016.
    """
    # Only the modules in PYTHON_PATH_MODULES use the Python path.
    module_name = request.module.__name__.rsplit(".", 1)[-1]
    if module_name in PYTHON_PATH_MODULES:
        # The native path does not call apply_to_points. Use Python here.
        monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
