"""Decide if LoopStructural uses the native fault kernel.

The native fault kernel is the Rust package ``loop_interpolation_core``.
"""

from __future__ import annotations

import os
from collections.abc import Sequence

import numpy as np

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


class _Unsupported(Exception):
    """Raised when a model part cannot go to the native kernel."""


def _grid(member: object, refs: list, grids: list) -> dict:
    """Return the grid dict of one frame member or unconformity field."""
    from loop_common.supports import StructuredGrid
    from loop_interpolation import DiscreteInterpolator

    if member.builder is not None:
        member.builder.up_to_date()
    interpolator = member.interpolator
    if not isinstance(interpolator, DiscreteInterpolator):
        raise _Unsupported("interpolator")
    for name in ("evaluate_value", "_evaluate_value_local", "_evaluate_gradient_local"):
        if getattr(type(interpolator), name) is not getattr(DiscreteInterpolator, name):
            raise _Unsupported("interpolator")
    support = interpolator.support
    if not isinstance(support, StructuredGrid) or not support._has_uniform_lookup():
        raise _Unsupported("support")
    # Same stale-frame check as the Python path.
    interpolator._check_frame()
    interpolator.update()
    box = interpolator.bounding_box
    grid = {
        "origin": np.asarray(support.origin, dtype=float).tolist(),
        "nsteps": [int(n) for n in support.nsteps],
        "step_vector": np.asarray(support.step_vector, dtype=float).tolist(),
        "local_origin": None if box is None else np.asarray(box.local_origin, dtype=float).tolist(),
        "local_rotation": None
        if box is None
        else np.asarray(box.local_rotation, dtype=float).tolist(),
        # The raw array is kept here; it is converted only when a kernel is built.
        "values": interpolator.c,
    }
    refs.append(interpolator.c)
    grids.append(grid)
    return grid


def _part(function: object) -> dict:
    """Return the profile dict of one profile function."""
    from .fault._fault_function import Composite, CubicFunction, Ones, Zeros, smooth_peak

    if function is smooth_peak:
        return {"kind": "smooth_peak"}
    kind = type(function)
    if kind is Ones:
        return {"kind": "ones"}
    if kind is Zeros:
        return {"kind": "zeros"}
    if kind is Composite:
        return {
            "kind": "composite",
            "positive": _part(function.positive),
            "negative": _part(function.negative),
        }
    if kind is CubicFunction:
        function.solve()
        return {
            "kind": "cubic",
            "w": [float(x) for x in function.w],
            "lim": [float(function.lim[0]), float(function.lim[1])],
            "min": float(function.min_v),
            "max": float(function.max_v),
            "value_points": [[float(x), float(y)] for x, y in function.value_points],
        }
    raise _Unsupported("profile")


class _Build:
    """State of one spec build: the chain faults, their specs and the arrays used."""

    def __init__(self) -> None:
        self.faults: list = []
        self.index: dict = {}
        self.specs: list = []
        self.refs: list = []
        self.grids: list = []

    def add(self, fault: object) -> int:
        """Return the chain index of a fault; append it if it is new."""
        i = self.index.get(id(fault))
        if i is None:
            i = len(self.faults)
            self.index[id(fault)] = i
            self.faults.append(fault)
        return i

    def indices(self, faults: object, owner: object) -> list:
        """Return the chain indices of the faults in restore order (cached on owner)."""
        from ...utils._fault_block import restore_order

        return [self.add(f) for f in restore_order(faults, owner=owner)]


def _region(region: object, build: _Build) -> dict:
    """Return the region dict of one region object."""
    from ...utils.regions import NegativeRegion, PositiveRegion
    from ._unconformity_feature import UnconformityFeature
    from .fault._fault_segment import FaultSegment

    if isinstance(region, UnconformityFeature):
        return {
            "kind": "unconformity",
            "field": _grid(region, build.refs, build.grids),
            "faults": build.indices(region.faults, region),
            "value": float(region.value),
            "sign": bool(region.sign),
        }
    if isinstance(region, (PositiveRegion, NegativeRegion)):
        if not isinstance(region.feature, FaultSegment) or region.vector is None:
            raise _Unsupported("region")
        if region.point is None:
            raise _Unsupported("region")
        return {
            "kind": "positive" if isinstance(region, PositiveRegion) else "negative",
            "fault": build.add(region.feature),
            "vector": [float(x) for x in region.vector],
            "point": [float(x) for x in region.point],
        }
    raise _Unsupported("region")


def _fault_spec(fault: object, build: _Build) -> dict:
    """Return the spec dict of one fault of the chain."""
    from ._geological_feature import GeologicalFeature
    from .fault._fault_function import FaultDisplacement
    from .fault._fault_segment import FaultSegment

    if not isinstance(fault, FaultSegment):
        raise _Unsupported("fault")
    if fault.builder is not None:
        fault.builder.up_to_date()
    function = fault.faultfunction
    if not isinstance(function, FaultDisplacement) or int(fault.steps) < 1:
        raise _Unsupported("profile")
    members = [fault[i] for i in range(3)]
    for member in members:
        if type(member).evaluate_value is not GeologicalFeature.evaluate_value:
            raise _Unsupported("member")
        if type(member).evaluate_gradient is not GeologicalFeature.evaluate_gradient:
            raise _Unsupported("member")
        same = len(member.faults) == len(members[0].faults) and all(
            a is b for a, b in zip(member.faults, members[0].faults)
        )
        if not same or (len(member.faults) > 0 and not member.faults_enabled):
            raise _Unsupported("member faults")
    return {
        "name": str(fault.name),
        "displacement": float(fault.displacement),
        "steps": int(fault.steps),
        "frame": [_grid(m, build.refs, build.grids) for m in members],
        "cutters": build.indices(members[0].faults, members[0]),
        "regions": [[_region(r, build) for r in m.regions] for m in members],
        "profile": {
            "scale": float(function.scale),
            "gx": _part(function.gx),
            "gy": _part(function.gy),
            "gz": _part(function.gz),
        },
    }


def _freeze(value: object) -> object:
    """Return a hashable copy of a spec part; an array of node values becomes its id."""
    if isinstance(value, dict):
        return tuple(
            (k, id(v) if k == "values" else _freeze(v)) for k, v in sorted(value.items())
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(v) for v in value)
    return value


def chain_for(feature: object) -> tuple | None:
    """Return the native fault kernel of a feature and its faults, or None.

    None means that the Python path must be used (FR-008). The kernel is kept
    in ``feature._native_chain`` and is made again only if the state token
    changes (FR-009, FR-010).

    Parameters
    ----------
    feature : BaseFeature
        The feature whose faults are restored.

    Returns
    -------
    tuple or None
        ``(FaultChain, faults, specs)`` where ``faults`` is the chain list and
        ``specs`` the spec dict of each fault, or None if any part is not
        supported.
    """
    from ...utils import _fault_block as fb

    if not native_enabled() or not getattr(feature, "faults_enabled", False):
        return None
    if not feature.faults:
        return None
    block = {
        "side_epsilon": fb.BLOCK_SIDE_EPSILON,
        "extension": fb.BLOCK_EXTENSION,
        "exact": fb.BLOCK_EXTENSION_EXACT,
        "limit": fb.BLOCK_EXTENSION_LIMIT,
        "step": fb.BLOCK_EXTENSION_STEP,
    }
    build = _Build()
    try:
        # Only the set of faults is needed here, so no order is computed.
        for f in fb.with_cutters(feature.faults):
            build.add(f)
        k = 0
        while k < len(build.faults):
            build.specs.append(_fault_spec(build.faults[k], build))
            k += 1
    except (_Unsupported, TypeError, ValueError, IndexError, KeyError, AttributeError):
        return None
    spec = {"version": 1, "block": block, "faults": build.specs}
    token = (tuple(block.items()), tuple(id(f) for f in build.faults), _freeze(build.specs))
    cached = getattr(feature, "_native_chain", None)
    if cached is not None and cached[0] == token:
        return cached[2], cached[3], cached[4]
    for grid in build.grids:
        grid["values"] = np.asarray(grid["values"], dtype=float)
    chain = _core.FaultChain(spec)
    try:
        feature._native_chain = (token, build.refs, chain, build.faults, build.specs)
    except AttributeError:
        pass
    return chain, build.faults, build.specs


def restore(
    feature: object,
    points: np.ndarray,
    *,
    reverse: bool = False,
    block: dict | None = None,
    restored_through: Sequence = (),
) -> np.ndarray | None:
    """Restore points through the faults of a feature with the native kernel.

    Parameters
    ----------
    feature : BaseFeature
        The feature whose faults are restored.
    points : np.ndarray
        Float array of shape (N, 3).
    reverse : bool
        Move the points forward instead of restoring them.
    block : dict, optional
        Fault block, see ``check_block``.
    restored_through : sequence of FaultSegment
        Faults whose restoration the points already have.

    Returns
    -------
    np.ndarray or None
        The restored points, or None if the Python path must be used.
    """
    from ...utils._fault_block import check_block, is_in, restore_order

    if not native_enabled():
        return None
    points = np.asarray(points)
    if points.dtype.kind != "f" or points.ndim != 2 or points.shape[1] != 3:
        return None
    found = chain_for(feature)
    if found is None:
        return None
    chain, faults, specs = found
    index = {id(f): i for i, f in enumerate(faults)}
    ordered = restore_order(feature.faults, reverse=reverse, owner=feature)
    order = [index[id(f)] for f in ordered if not is_in(f, restored_through)]
    named = check_block(block) or {}
    sides = [named.get(f.name, 2) for f in faults]
    restored = [is_in(f, restored_through) for f in faults]
    # A cutter that is not restored at the turn of its fault needs the Python path.
    done = {i for i, flag in enumerate(restored) if flag}
    for i in order:
        if not all(c in done for c in specs[i]["cutters"]):
            return None
        done.add(i)
    points = np.ascontiguousarray(points, dtype=float)
    return chain.restore(
        points,
        np.asarray(order, dtype=np.int64),
        bool(reverse),
        np.asarray(sides, dtype=np.int8),
        np.asarray(restored, dtype=bool),
        native_threads(),
    )
