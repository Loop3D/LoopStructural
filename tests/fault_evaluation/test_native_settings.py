"""Tests for the native fault kernel settings (FR-006, FR-006a, FR-006b, FR-012).

The tests use a fake ``loop_model_core`` module. They check the setting,
the environment variable, the version check and the thread count of
``LoopStructural.modelling.features._native_faults``.

The module is reloaded before each test. This resets its one-time warning flag.
"""

from __future__ import annotations

import importlib
import logging
import sys
from types import SimpleNamespace

import pytest

from LoopStructural import LoopStructuralConfig

MODULE_NAME = "LoopStructural.modelling.features._native_faults"
KERNEL_NAME = "loop_model_core"
ENV_VAR = "LOOPSTRUCTURAL_NATIVE_FAULTS"


def _set_kernel_version(monkeypatch, version):
    """Install a fake kernel module with the given ``FAULT_KERNEL_VERSION``.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture used to undo the change after the test.
    version : int
        The version number the fake kernel reports.
    """
    monkeypatch.setitem(sys.modules, KERNEL_NAME, SimpleNamespace(FAULT_KERNEL_VERSION=version))


@pytest.fixture
def native(monkeypatch):
    """Return a freshly reloaded ``_native_faults`` module with clean settings.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture used to set and undo the settings, the environment and ``sys.modules``.

    Returns
    -------
    module
        The reloaded ``_native_faults`` module.
    """
    # The settings do not exist yet, so raising=False is needed.
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_threads", None, raising=False)
    monkeypatch.delenv(ENV_VAR, raising=False)
    _set_kernel_version(monkeypatch, 1)
    module = importlib.import_module(MODULE_NAME)
    importlib.reload(module)
    yield module
    # Drop the fake kernel from the module state, so other tests in the same worker do not use it.
    importlib.reload(module)


@pytest.fixture
def module_log(caplog):
    """Capture the warnings that the ``_native_faults`` logger emits.

    The logger does not always reach the root logger. The test therefore
    attaches the caplog handler to the logger itself.

    Parameters
    ----------
    caplog : pytest.LogCaptureFixture
        Pytest fixture that collects the log records.

    Yields
    ------
    logging.Logger
        The logger of the module under test.
    """
    logger = logging.getLogger(MODULE_NAME)
    old_propagate = logger.propagate
    logger.propagate = False
    logger.addHandler(caplog.handler)
    caplog.set_level(logging.WARNING, logger=MODULE_NAME)
    try:
        yield logger
    finally:
        logger.removeHandler(caplog.handler)
        logger.propagate = old_propagate


@pytest.mark.parametrize("value", ["0", "false", "OFF", "no"])
def test_env_variable_disables_native(monkeypatch, native, value):
    """The environment variable turns the native kernel off for these values."""
    monkeypatch.setenv(ENV_VAR, value)
    assert native.native_enabled() is False


def test_unset_variable_enables_native(native):
    """With no setting and no variable, the native kernel is on."""
    assert native.native_enabled() is True


def test_setting_true_wins_over_env_zero(monkeypatch, native):
    """A True setting wins over the environment variable value "0"."""
    monkeypatch.setenv(ENV_VAR, "0")
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", True, raising=False)
    assert native.native_enabled() is True


def test_setting_false_disables_native(monkeypatch, native):
    """A False setting turns the native kernel off."""
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False, raising=False)
    assert native.native_enabled() is False


def test_import_failure_disables_native_without_warning(monkeypatch, native, module_log, caplog):
    """A kernel that cannot be imported gives False and no warning."""
    monkeypatch.setitem(sys.modules, KERNEL_NAME, None)
    assert native.native_enabled() is False
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert warnings == []


def test_version_mismatch_disables_native_and_warns_once(monkeypatch, native, module_log, caplog):
    """A version 2 kernel gives False, and the warning is logged only once."""
    _set_kernel_version(monkeypatch, 2)
    assert native.native_enabled() is False
    assert native.native_enabled() is False
    matches = [
        r
        for r in caplog.records
        if r.levelno == logging.WARNING and "2" in r.getMessage() and "1" in r.getMessage()
    ]
    assert len(matches) == 1


def test_native_threads_none_is_zero(native):
    """With no thread setting, native_threads() returns 0."""
    assert native.native_threads() == 0


def test_native_threads_one(monkeypatch, native):
    """A thread setting of 1 is returned as 1."""
    monkeypatch.setattr(LoopStructuralConfig, "native_threads", 1, raising=False)
    assert native.native_threads() == 1


def test_negative_thread_count_means_all_cores(monkeypatch, native):
    """A negative thread setting means all cores, as -1 in joblib, so it is returned as 0."""
    monkeypatch.setattr(LoopStructuralConfig, "native_threads", -1, raising=False)
    assert native.native_threads() == 0


def test_subclass_of_fault_displacement_uses_python_path(monkeypatch):
    """A subclass of FaultDisplacement takes the Python path (FR-008)."""
    pytest.importorskip("loop_model_core")
    from LoopStructural.modelling.features import _native_faults
    from LoopStructural.modelling.features.fault._fault_function import FaultDisplacement
    from tests.fault_evaluation.models import build_model

    # The fake kernel of the native fixture is undone after its test, so the real kernel is used.
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", True, raising=False)
    monkeypatch.setattr(LoopStructuralConfig, "native_threads", None, raising=False)
    monkeypatch.delenv(ENV_VAR, raising=False)
    importlib.reload(_native_faults)

    class Sub(FaultDisplacement):
        pass

    model = build_model("f1_u0")
    feature = next(f for f in model.features if getattr(f, "faults", None))
    fault = feature.faults[0]
    old = fault.faultfunction
    sub = Sub.__new__(Sub)
    sub.__dict__.update(old.__dict__)
    fault.faultfunction = sub
    assert fault.faultfunction is sub
    assert _native_faults.chain_for(feature) is None
