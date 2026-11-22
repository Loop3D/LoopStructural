"""Adapters that turn outside objects into loop_link scene objects."""

from __future__ import annotations

from typing import Any, get_args

from loop_link.errors import InvalidInputError
from loop_link.objects import SceneObject

_ACCEPTED_INPUTS = (
    "Accepted inputs: loop_link objects (Points, Vectors, Lines, Surface, BlockModel, "
    "TetraMesh, Drillhole); pyvista meshes and grids; LoopStructural Surface and "
    "StructuredGrid; pandas tables for drillholes "
    "(use loop_link.adapters.pandas.drillholes_from_pandas)."
)


def convert(obj: Any, **hints: Any) -> SceneObject:  # noqa: ANN401
    """Convert an input object to a loop_link scene object.

    Parameters
    ----------
    obj : Any
        A loop_link object, a pyvista mesh or grid, or a LoopStructural object.
    **hints : Any
        Options for the pyvista adapter.

    Returns
    -------
    SceneObject
        ``obj`` itself when it is a scene object, else the converted object.

    Raises
    ------
    loop_link.errors.InvalidInputError
        The input type is not accepted.
    ImportError
        The adapter module for the input type is not available.
    """
    if isinstance(obj, get_args(SceneObject)):
        return obj
    package = type(obj).__module__.split(".", 1)[0]
    if package == "pyvista":
        from loop_link.adapters.pyvista import from_pyvista

        return from_pyvista(obj, **hints)
    if package == "loop_common":
        from loop_link.adapters.loopstructural import from_loopstructural

        return from_loopstructural(obj)
    raise InvalidInputError(_ACCEPTED_INPUTS)
