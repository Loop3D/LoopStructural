"""Decide if LoopStructural uses the native fault kernel.

The native fault kernel is the Rust package ``loop_interpolation_core``.
"""

from __future__ import annotations

import os

from ...utils import getLogger

logger = getLogger(__name__)

REQUIRED_FAULT_KERNEL_VERSION: int = 1

_core = None
_core_checked = False
_version_warned = False


def native_enabled() -> bool:
    """Return True if the native fault kernel can be used and is not turned off.

    The setting ``LoopStructuralConfig.native_faults`` has priority over the
    environment variable ``LOOPSTRUCTURAL_NATIVE_FAULTS``. The kernel import
    and the version check run once for each session.

    Returns
    -------
    bool
        True if the kernel is on, installed and has the required version.
    """
    global _core, _core_checked, _version_warned
    from LoopStructural import LoopStructuralConfig

    setting = LoopStructuralConfig.native_faults
    if setting is None:
        if os.environ.get("LOOPSTRUCTURAL_NATIVE_FAULTS", "").lower() in {
            "0",
            "false",
            "off",
            "no",
        }:
            return False
    elif not setting:
        return False
    if not _core_checked:
        _core_checked = True
        try:
            import loop_interpolation_core

            _core = loop_interpolation_core
        except ImportError:
            _core = None
    if _core is None:
        return False
    found = getattr(_core, "FAULT_KERNEL_VERSION", None)
    if found != REQUIRED_FAULT_KERNEL_VERSION:
        if not _version_warned:
            _version_warned = True
            logger.warning(
                f"native fault kernel version {found} found, LoopStructural needs "
                f"version {REQUIRED_FAULT_KERNEL_VERSION}; using the Python path"
            )
        return False
    return True


def native_threads() -> int:
    """Return the thread count for the native fault kernel.

    Returns
    -------
    int
        The setting ``LoopStructuralConfig.native_threads``, or 0 if it is None.
    """
    from LoopStructural import LoopStructuralConfig

    return LoopStructuralConfig.native_threads or 0
