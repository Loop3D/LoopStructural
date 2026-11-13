"""Helpers for evaluating a field in one fault block.

A fault block gives, for some faults, the side of the fault that every point
is taken to be on:

* ``+1``: the hanging wall (the side where the fault frame coordinate gx > 0),
* ``-1``: the footwall (gx < 0),
* ``0``: the fault is absent: it moves no point.

A fault that the block does not name keeps the side of each point. The block
is passed down the evaluation as an argument. No object is changed, so two
evaluations with different blocks can run at the same time.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

BLOCK_SIDE_EPSILON = 1e-9
"""Smallest absolute value of a forced fault frame coordinate.

A point exactly on the fault (gx = 0) is moved to this distance on the
forced side, so that the displacement profile takes the value of that side.
"""


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
