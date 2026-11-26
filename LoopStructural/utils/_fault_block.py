"""Helpers for evaluating a field in one fault block.

A fault block gives, for some faults, the side of the fault that every point
is taken to be on:

* ``+1``: the hanging wall (the side where the fault frame coordinate gx > 0),
* ``-1``: the footwall (gx < 0),
* ``0``: the fault is absent: it moves no point.

A fault that the block does not name keeps the side of each point. The block
is passed down the evaluation as an argument. No object is changed, so two
evaluations with different blocks can run at the same time.

**The block on the other side of its fault.** A block field is also evaluated
at points on the other side of the fault. There the displacement of the fault
is the smooth continuation of the displacement profile of the block's own
side across gx = 0 (``BLOCK_EXTENSION = "smooth"``, see
``LoopStructural.modelling.features.fault._fault_function.block_displacement``):
the value, the slope and the curvature of the profile at the fault are kept,
so the block field is C2 across the fault surface and has no sharp change of
gradient there. The continuation has no geological meaning; a mesher needs
it only within a few cells of the fault.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np

BLOCK_SIDE_EPSILON = 1e-9
"""Smallest absolute value of a forced fault frame coordinate.

A point exactly on the fault (gx = 0) is moved to this distance on the
forced side, so that the displacement profile takes the value of that side.
"""


BLOCK_EXTENSION = "smooth"
"""How the displacement of a block continues on the other side of its fault.

``"smooth"`` (default): the continuation of the profile of the block's side.
``"mirror"``: the mirror image, gx becomes ``side * |gx|`` (the first draft
and the fault benchmark study data). The mirror image has a kink where the
profile has a slope at the fault. It is kept only to compare with data made
with it. This is a module setting, read at each evaluation: do not change it
while other threads evaluate.
"""

BLOCK_EXTENSION_EXACT = 0.25
"""Distance from the fault, in units of gx, up to which the continuation is
the exact continuation of the profile."""

BLOCK_EXTENSION_LIMIT = 0.5
"""Distance from the fault, in units of gx, beyond which the continuation is
constant. Between the two distances it is blended (C2), so the continuation
cannot grow without limit. The library profiles decay over |gx| = 1."""

BLOCK_EXTENSION_STEP = 1e-3
"""Step in gx of the one-sided finite differences that give the slope and
the curvature at the fault of a profile that is not a polynomial."""


def check_block(block) -> dict[str, int] | None:
    """Check a fault block and return it as a dictionary of fault names.

    Parameters
    ----------
    block : Mapping or None
        Maps a fault (its name or the fault object) to +1, -1 or 0.

    Returns
    -------
    dict of str to int, or None
        The block with fault names as keys, or None if ``block`` is None.

    Raises
    ------
    TypeError
        If ``block`` is not a mapping.
    ValueError
        If a side is not +1, -1 or 0.
    """
    if block is None:
        return None
    if not isinstance(block, Mapping):
        raise TypeError("block must map fault names to +1, -1 or 0")
    out = {}
    for fault, side in block.items():
        name = fault if isinstance(fault, str) else getattr(fault, "name", None)
        if not isinstance(name, str):
            raise TypeError(f"block key {fault!r} is not a fault or a fault name")
        if isinstance(side, bool) or side not in (-1, 0, 1):
            raise ValueError(f"side of fault {name} must be +1, -1 or 0, not {side!r}")
        out[name] = int(side)
    return out


def block_side(block: Mapping[str, int] | None, name: str) -> int | None:
    """Return the side of one fault in a block, or None if the block does not name it."""
    if not block:
        return None
    return block.get(name)


def force_side(values: np.ndarray, side: int | None) -> np.ndarray:
    """Put fault frame values on one side of the fault.

    The value v becomes ``side * max(|v|, BLOCK_SIDE_EPSILON)``: the mirror
    image of the field on the other side. The result is continuous across the
    fault. NaN values (no support, or outside a region of the fault) stay NaN.

    This decides the side of a point (regions, abutting faults). The
    displacement of a block on the other side of its fault does not use the
    mirror image, see ``BLOCK_EXTENSION``.

    Parameters
    ----------
    values : np.ndarray
        Fault frame coordinate gx.
    side : int or None
        +1 or -1. Any other value returns the values unchanged.

    Returns
    -------
    np.ndarray
        A new array.
    """
    out = np.array(values, dtype=float, copy=True)
    if side not in (-1, 1):
        return out
    known = ~np.isnan(out)
    out[known] = side * np.maximum(np.abs(out[known]), BLOCK_SIDE_EPSILON)
    return out


def block_kwargs(block: Mapping[str, int] | None) -> dict:
    """Keyword arguments that pass a block on, or none if there is no block.

    Callers use this so that an evaluation without a block calls the same
    signatures as before.
    """
    return {} if block is None else {"block": block}


def is_in(fault, faults) -> bool:
    """True if ``fault`` is one of ``faults`` (by identity)."""
    return any(fault is f for f in faults)


def _cutters(fault: object) -> list:
    """Return the faults that displace the frame of ``fault`` (its ``faults=[...]``).

    Before the frame is built only its builders know them; after the build the
    frame members do.
    """
    found = []
    builder = getattr(fault, "builder", None)
    sources = [getattr(b, "faults", None) or [] for b in getattr(builder, "builders", None) or []]
    try:
        sources += [getattr(fault[i], "faults", None) or [] for i in range(3)]
    except (IndexError, KeyError, TypeError):
        pass
    for source in sources:
        for cutter in source:
            if not is_in(cutter, found):
                found.append(cutter)
    return found


def cutters_first(faults: Iterable) -> list:
    """Order ``faults`` so that a fault comes after every fault of the list that cuts it.

    A fault created with ``faults=[A]`` is faulted by the younger fault ``A``, so
    ``A`` must be restored first. Faults without such a relation keep the order
    of the list. A cycle of relations keeps the order of the list for the faults
    in it.
    """
    faults = list(faults)
    pending, ordered = list(faults), []
    while pending:
        for fault in pending:
            if not any(is_in(c, pending) and c is not fault for c in _cutters(fault)):
                break
        else:
            fault = pending[0]  # a cycle: keep the order of the list
        ordered.append(fault)
        pending = [f for f in pending if f is not fault]
    return ordered


def restore_through(points, faults, *, reverse=False, block=None, restored_through=()):
    """Restore points through a list of faults, each fault once.

    The faults are applied in the order of the list (a feature keeps its
    faults youngest first), except that a fault faulted by another fault of the
    list (``faults=[...]``) comes after it; when the points move forward
    (``reverse``) the cutting fault comes last. Each fault is told which faults the points have
    already been restored through, so that its frame, which was built in the
    space restored by the younger faults that displace it, is evaluated at the
    points as they are and does not apply those younger faults a second time.

    Parameters
    ----------
    points : np.ndarray
        Nx3 locations
    faults : list of FaultSegment
        faults in the order to apply them
    reverse : bool
        move the points forward instead of restoring them
    block : dict, optional
        checked fault block
    restored_through : sequence of FaultSegment
        faults whose restoration the points already have; they are skipped

    Returns
    -------
    np.ndarray
        the restored points
    """
    done = list(restored_through)
    if reverse:
        # Forward movement: a cutting fault moves after the faults it cuts. Without
        # such relations the order of the list does not change.
        ordered = cutters_first(list(faults)[::-1])[::-1]
    else:
        ordered = cutters_first(faults)
    for fault in ordered:
        if is_in(fault, done):
            continue
        # Pass only the keywords that are needed, so a fault-like object with
        # the plain apply_to_points(points) signature still works.
        kwargs = block_kwargs(block)
        if reverse:
            kwargs["reverse"] = True
        if done:
            kwargs["restored_through"] = tuple(done)
        points = fault.apply_to_points(points, **kwargs)
        done.append(fault)
    return points
