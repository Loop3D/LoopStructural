"""Adapter that turns LoopStructural geometry objects into loop_link scene objects."""

from __future__ import annotations

import sys
from typing import Any

import numpy as np

from loop_link.errors import InvalidInputError, MissingExtraError
from loop_link.objects import BlockModel, SceneObject, Surface


def from_loopstructural(obj: Any) -> SceneObject:  # noqa: ANN401
    """Convert a LoopStructural Surface or StructuredGrid to a scene object.

    The coordinates stay in world coordinates. No other transform is applied.

    Parameters
    ----------
    obj : Any
        A ``loop_common.geometry.Surface`` or ``loop_common.geometry.StructuredGrid``.

    Returns
    -------
    SceneObject
        A loop_link ``Surface`` for a Surface, or a ``BlockModel`` for a StructuredGrid.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The input is not a Surface or a StructuredGrid.
    loop_link.errors.MissingExtraError
        The ``loop_common`` package is not installed.
    """
    from loop_link.adapters import _ACCEPTED_INPUTS

    module = sys.modules.get("loop_common.geometry")
    if module is None:
        try:
            import loop_common.geometry as module
        except ImportError as error:
            raise MissingExtraError("loopstructural") from error
    if isinstance(obj, module.Surface):
        return Surface(
            vertices=obj.vertices,
            triangles=obj.triangles,
            properties=dict(obj.properties or {}),
            cell_properties=dict(obj.cell_properties or {}),
        )
    if isinstance(obj, module.StructuredGrid):
        return BlockModel(
            origin=obj.origin,
            cell_sizes=obj.step_vector,
            cells=obj.nsteps_cells,
            properties=_flatten(obj.properties),
            cell_properties=_flatten(obj.cell_properties),
        )
    raise InvalidInputError(_ACCEPTED_INPUTS)


def _flatten(properties: dict[str, Any] | None) -> dict[str, Any]:
    """Return the grid properties as 1-D arrays in x-fastest order.

    ``StructuredGrid.vtk()`` reshapes each property to ``(n, -1)`` with ``order="F"``,
    so a property can be a grid of shape ``(nx, ny, nz)`` or have more than one value
    for each node, for example a vector of shape ``(n, 3)``. A grid is flattened in
    x-fastest order. Each column of a 2-D array becomes one property with the name
    ``"<name>_<column>"``. A 2-D array with one column keeps its name.
    """
    result: dict[str, Any] = {}
    for name, values in (properties or {}).items():
        array = np.asarray(values)
        if array.ndim == 3:
            array = array.ravel(order="F")
        if array.ndim == 2 and array.shape[1] == 1:
            array = array[:, 0]
        if array.ndim == 2:
            for column in range(array.shape[1]):
                result[f"{name}_{column}"] = array[:, column]
            continue
        result[name] = array
    return result
