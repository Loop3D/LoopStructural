"""The object model of loop_link.

Each class is a frozen dataclass. A class checks its input in ``__post_init__``
and stores the arrays as numpy arrays. A check that fails raises
:class:`loop_link.errors.InvalidInputError` with the name of the field.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Union

import numpy as np

from loop_link.errors import InvalidInputError
from loop_link.orientation import directions_to_lineation, normals_to_orientation

__all__ = [
    "DEFAULT_CATEGORY_COLOURS",
    "ArrayLike",
    "BlockModel",
    "Categorical",
    "Drillhole",
    "Lines",
    "Points",
    "Property",
    "SceneObject",
    "Surface",
    "TetraMesh",
    "Vectors",
]

ArrayLike = Union[np.ndarray, Sequence]

#: The ten tab10 colours as RGBA, with alpha 255.
DEFAULT_CATEGORY_COLOURS: tuple[tuple[int, int, int, int], ...] = (
    (31, 119, 180, 255),
    (255, 127, 14, 255),
    (44, 160, 44, 255),
    (214, 39, 40, 255),
    (148, 103, 189, 255),
    (140, 86, 75, 255),
    (227, 119, 194, 255),
    (127, 127, 127, 255),
    (188, 189, 34, 255),
    (23, 190, 207, 255),
)

_MAX_CODES = 2**32


def _vertices(value: object, name: str = "vertices") -> np.ndarray:
    """Return vertices as a finite float64 array with shape (n, 3)."""
    try:
        arr = np.array(value, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise InvalidInputError(f"The field '{name}' must be numbers. Check the input.") from exc
    if arr.size == 0 and arr.ndim == 1:
        arr = arr.reshape(0, 3)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise InvalidInputError(f"The field '{name}' must have shape (n, 3). Check the input.")
    if not np.all(np.isfinite(arr)):
        raise InvalidInputError(f"The field '{name}' must be finite numbers. Check the input.")
    return arr


def _vector3(value: object, name: str) -> np.ndarray:
    """Return a finite float64 array with shape (3,)."""
    try:
        arr = np.array(value, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise InvalidInputError(f"The field '{name}' must be numbers. Check the input.") from exc
    if arr.shape != (3,):
        raise InvalidInputError(f"The field '{name}' must have shape (3,). Check the input.")
    if not np.all(np.isfinite(arr)):
        raise InvalidInputError(f"The field '{name}' must be finite numbers. Check the input.")
    return arr


def _integer_array(value: object, name: str) -> np.ndarray:
    """Return ``value`` as an int64 array.

    Raises
    ------
    InvalidInputError
        If the values are not integers.
    """
    try:
        arr = np.array(value, copy=True)
    except (TypeError, ValueError) as exc:
        raise InvalidInputError(f"The field '{name}' must be integers. Check the input.") from exc
    if arr.size == 0:
        return arr.astype(np.int64)
    if arr.dtype.kind not in "iu":
        raise InvalidInputError(f"The field '{name}' must be integers. Check the input.")
    return arr.astype(np.int64)


def _index_array(value: object, name: str, width: int, count: int) -> np.ndarray:
    """Return an index array with shape (m, width) and values in [0, count)."""
    arr = _integer_array(value, name)
    if arr.size == 0 and arr.ndim == 1:
        arr = arr.reshape(0, width)
    if arr.ndim != 2 or arr.shape[1] != width:
        raise InvalidInputError(
            f"The field '{name}' must have shape (m, {width}). Check the input."
        )
    if arr.size and (arr.min() < 0 or arr.max() >= count):
        raise InvalidInputError(
            f"The field '{name}' has an index outside the range 0 to {count - 1}. "
            "Check the input."
        )
    return arr


@dataclass(frozen=True, eq=False)
class Categorical:
    """A property that holds categories.

    Parameters
    ----------
    codes : array-like of int
        One code for each item. A code is an index into ``names``.
    names : sequence of str
        The name of each category.
    colours : sequence of tuple, optional
        One RGBA colour for each category. A category without a colour uses
        ``DEFAULT_CATEGORY_COLOURS[i % 10]``.
    """

    codes: ArrayLike
    names: Sequence[str]
    colours: Sequence[tuple[int, int, int, int]] | None = None

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        codes = _integer_array(self.codes, "codes")
        if codes.ndim != 1:
            raise InvalidInputError("The field 'codes' must be one-dimensional. Check the input.")
        limit = min(len(self.names), _MAX_CODES)
        if codes.size and (codes.min() < 0 or codes.max() >= limit):
            raise InvalidInputError(
                f"The field 'codes' must have values from 0 to {limit - 1}. Check the input."
            )
        colours = self.colours
        if colours is not None:
            try:
                arr = np.array(colours, copy=True)
            except (TypeError, ValueError) as exc:
                raise InvalidInputError(
                    "The field 'colours' must be RGBA numbers. Check the input."
                ) from exc
            if arr.shape != (len(self.names), 4) or arr.dtype.kind not in "iu":
                raise InvalidInputError(
                    "The field 'colours' must have one RGBA tuple of 4 integers for each name. "
                    "Check the input."
                )
            if arr.size and (arr.min() < 0 or arr.max() > 255):
                raise InvalidInputError(
                    "The field 'colours' must have values from 0 to 255. Check the input."
                )
            colours = tuple(tuple(int(v) for v in row) for row in arr)
        object.__setattr__(self, "codes", codes)
        object.__setattr__(self, "colours", colours)


Property = Union[np.ndarray, Sequence, Categorical]


def _properties(properties: dict[str, Any] | None, length: int, field_name: str) -> dict[str, Any]:
    """Check a mapping of properties and return a new mapping.

    Parameters
    ----------
    properties : dict
        Property name to values or :class:`Categorical`.
    length : int
        The required number of values for each property.
    field_name : str
        The name of the field. The name is used in error messages.

    Returns
    -------
    dict
        The properties as numpy arrays. A :class:`Categorical` stays the same object.

    Raises
    ------
    InvalidInputError
        If a property is not one-dimensional or has the wrong length.
    """
    if properties is None:
        properties = {}
    if not isinstance(properties, Mapping):
        raise InvalidInputError(f"The field '{field_name}' must be a mapping. Check the input.")
    result: dict[str, Any] = {}
    for key, value in properties.items():
        if isinstance(value, Categorical):
            size = len(value.codes)
            result[key] = value
        else:
            try:
                arr = np.array(value, copy=True)
            except (TypeError, ValueError) as exc:
                raise InvalidInputError(
                    f"The property '{key}' in '{field_name}' is not valid. Check the input."
                ) from exc
            if arr.ndim != 1:
                raise InvalidInputError(
                    f"The property '{key}' in '{field_name}' must be one-dimensional. "
                    "Check the input."
                )
            size = arr.shape[0]
            result[key] = arr
        if size != length:
            raise InvalidInputError(
                f"The property '{key}' in '{field_name}' has {size} values, but {length} "
                "values are needed. Check the input."
            )
    return result


@dataclass(frozen=True, eq=False)
class Points:
    """A set of points.

    Parameters
    ----------
    vertices : array-like
        Coordinates with shape (n, 3).
    properties : dict
        Property name to n values or a :class:`Categorical`.
    """

    vertices: ArrayLike
    properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        verts = _vertices(self.vertices)
        object.__setattr__(self, "vertices", verts)
        object.__setattr__(
            self, "properties", _properties(self.properties, len(verts), "properties")
        )


@dataclass(frozen=True, eq=False)
class Vectors:
    """A set of points with orientation or lineation data.

    Parameters
    ----------
    vertices : array-like
        Coordinates with shape (n, 3).
    orientation : tuple of array-like, optional
        ``(dip_direction, dip)`` in degrees, n values each.
    lineation : tuple of array-like, optional
        ``(trend, plunge)`` in degrees, n values each.
    properties : dict
        Property name to n values or a :class:`Categorical`.
    """

    vertices: ArrayLike
    orientation: tuple[ArrayLike, ArrayLike] | None = None
    lineation: tuple[ArrayLike, ArrayLike] | None = None
    properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        verts = _vertices(self.vertices)
        object.__setattr__(self, "vertices", verts)
        for name in ("orientation", "lineation"):
            pair = getattr(self, name)
            if pair is None:
                continue
            try:
                first, second = pair
            except (TypeError, ValueError) as exc:
                raise InvalidInputError(
                    f"The field '{name}' must be a tuple of two arrays. Check the input."
                ) from exc
            arrays = []
            for part in (first, second):
                try:
                    arr = np.array(part, dtype=np.float64, copy=True)
                except (TypeError, ValueError) as exc:
                    raise InvalidInputError(
                        f"The field '{name}' must be numbers. Check the input."
                    ) from exc
                if arr.shape != (len(verts),) or not np.all(np.isfinite(arr)):
                    raise InvalidInputError(
                        f"The field '{name}' must have {len(verts)} finite values in each array. "
                        "Check the input."
                    )
                arrays.append(arr)
            object.__setattr__(self, name, (arrays[0], arrays[1]))
        object.__setattr__(
            self, "properties", _properties(self.properties, len(verts), "properties")
        )

    @classmethod
    def from_directions(
        cls,
        vertices: ArrayLike,
        directions: ArrayLike,
        *,
        kind: str,
        properties: dict[str, Property] | None = None,
    ) -> Vectors:
        """Make vectors from direction vectors.

        Parameters
        ----------
        vertices : array-like
            Coordinates with shape (n, 3).
        directions : array-like
            Vectors with shape (n, 3). They are plane normals or line directions.
        kind : {"normal", "lineation"}
            The meaning of ``directions``.
        properties : dict, optional
            Property name to n values or a :class:`Categorical`.

        Returns
        -------
        Vectors
            The vectors with ``orientation`` or ``lineation`` set.

        Raises
        ------
        InvalidInputError
            If ``kind`` is not "normal" or "lineation".
        """
        props = properties if properties is not None else {}
        if kind == "normal":
            return cls(vertices, orientation=normals_to_orientation(directions), properties=props)
        if kind == "lineation":
            return cls(vertices, lineation=directions_to_lineation(directions), properties=props)
        raise InvalidInputError("The argument 'kind' must be 'normal' or 'lineation'.")


@dataclass(frozen=True, eq=False)
class Lines:
    """A set of line segments.

    Parameters
    ----------
    vertices : array-like
        Coordinates with shape (n, 3).
    segments : array-like of int
        Vertex indices with shape (m, 2).
    properties : dict
        Vertex properties, n values each.
    cell_properties : dict
        Segment properties, m values each.
    """

    vertices: ArrayLike
    segments: ArrayLike
    properties: dict[str, Property] = field(default_factory=dict)
    cell_properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        verts = _vertices(self.vertices)
        segs = _index_array(self.segments, "segments", 2, len(verts))
        object.__setattr__(self, "vertices", verts)
        object.__setattr__(self, "segments", segs)
        object.__setattr__(
            self, "properties", _properties(self.properties, len(verts), "properties")
        )
        object.__setattr__(
            self,
            "cell_properties",
            _properties(self.cell_properties, len(segs), "cell_properties"),
        )


@dataclass(frozen=True, eq=False)
class Surface:
    """A triangle mesh.

    Parameters
    ----------
    vertices : array-like
        Coordinates with shape (n, 3).
    triangles : array-like of int
        Vertex indices with shape (m, 3).
    properties : dict
        Vertex properties, n values each.
    cell_properties : dict
        Triangle properties, m values each.
    """

    vertices: ArrayLike
    triangles: ArrayLike
    properties: dict[str, Property] = field(default_factory=dict)
    cell_properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        verts = _vertices(self.vertices)
        tris = _index_array(self.triangles, "triangles", 3, len(verts))
        object.__setattr__(self, "vertices", verts)
        object.__setattr__(self, "triangles", tris)
        object.__setattr__(
            self, "properties", _properties(self.properties, len(verts), "properties")
        )
        object.__setattr__(
            self,
            "cell_properties",
            _properties(self.cell_properties, len(tris), "cell_properties"),
        )


@dataclass(frozen=True, eq=False)
class BlockModel:
    """A regular grid of cells.

    Node and cell values are in x-fastest order.

    Parameters
    ----------
    origin : array-like
        The corner of the grid, shape (3,).
    cell_sizes : array-like
        The size of a cell on each axis, shape (3,). All values are above zero.
    cells : array-like of int
        The number of cells on each axis, shape (3,). All values are at least 1.
    properties : dict
        Node properties, (nx+1)*(ny+1)*(nz+1) values each.
    cell_properties : dict
        Cell properties, nx*ny*nz values each.
    """

    origin: ArrayLike
    cell_sizes: ArrayLike
    cells: ArrayLike
    properties: dict[str, Property] = field(default_factory=dict)
    cell_properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        origin = _vector3(self.origin, "origin")
        sizes = _vector3(self.cell_sizes, "cell_sizes")
        if np.any(sizes <= 0):
            raise InvalidInputError("The field 'cell_sizes' must be above zero. Check the input.")
        try:
            raw = np.asarray(self.cells)
        except (TypeError, ValueError) as exc:
            raise InvalidInputError("The field 'cells' must be integers. Check the input.") from exc
        if raw.shape != (3,):
            raise InvalidInputError("The field 'cells' must have shape (3,). Check the input.")
        if raw.dtype.kind == "f":
            if not np.all(np.isfinite(raw)) or np.any(raw != np.round(raw)):
                raise InvalidInputError("The field 'cells' must be integers. Check the input.")
            raw = raw.astype(np.int64)
        cells = _integer_array(raw, "cells")
        if np.any(cells < 1):
            raise InvalidInputError("The field 'cells' must be at least 1. Check the input.")
        nodes = int(np.prod(cells + 1))
        count = int(np.prod(cells))
        object.__setattr__(self, "origin", origin)
        object.__setattr__(self, "cell_sizes", sizes)
        object.__setattr__(self, "cells", cells)
        object.__setattr__(self, "properties", _properties(self.properties, nodes, "properties"))
        object.__setattr__(
            self,
            "cell_properties",
            _properties(self.cell_properties, count, "cell_properties"),
        )


@dataclass(frozen=True, eq=False)
class TetraMesh:
    """A mesh of tetrahedra.

    Parameters
    ----------
    vertices : array-like
        Coordinates with shape (n, 3).
    tetrahedra : array-like of int
        Vertex indices with shape (m, 4).
    properties : dict
        Vertex properties, n values each.
    cell_properties : dict
        Tetrahedron properties, m values each.
    """

    vertices: ArrayLike
    tetrahedra: ArrayLike
    properties: dict[str, Property] = field(default_factory=dict)
    cell_properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        verts = _vertices(self.vertices)
        tets = _index_array(self.tetrahedra, "tetrahedra", 4, len(verts))
        object.__setattr__(self, "vertices", verts)
        object.__setattr__(self, "tetrahedra", tets)
        object.__setattr__(
            self, "properties", _properties(self.properties, len(verts), "properties")
        )
        object.__setattr__(
            self,
            "cell_properties",
            _properties(self.cell_properties, len(tets), "cell_properties"),
        )


@dataclass(frozen=True, eq=False)
class Drillhole:
    """A drillhole with a survey and intervals.

    Parameters
    ----------
    collar : array-like
        The collar position, shape (3,).
    survey : array-like
        Survey stations with shape (k, 3).
    intervals : array-like
        Interval start and end with shape (j, 2).
    end_of_hole : float, optional
        The depth of the end of the hole.
    properties : dict
        Survey station properties, k values each.
    interval_properties : dict
        Interval properties, j values each.
    """

    collar: ArrayLike
    survey: ArrayLike
    intervals: ArrayLike
    end_of_hole: float | None = None
    properties: dict[str, Property] = field(default_factory=dict)
    interval_properties: dict[str, Property] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Check the input and store it as numpy arrays."""
        collar = _vector3(self.collar, "collar")
        survey = _vertices(self.survey, "survey")
        try:
            intervals = np.array(self.intervals, dtype=np.float64, copy=True)
        except (TypeError, ValueError) as exc:
            raise InvalidInputError(
                "The field 'intervals' must be numbers. Check the input."
            ) from exc
        if intervals.ndim != 2 or intervals.shape[1] != 2:
            raise InvalidInputError(
                "The field 'intervals' must have shape (j, 2). Check the input."
            )
        if not np.all(np.isfinite(intervals)):
            raise InvalidInputError(
                "The field 'intervals' must be finite numbers. Check the input."
            )
        if self.end_of_hole is not None:
            try:
                eoh = float(self.end_of_hole)
            except (TypeError, ValueError) as exc:
                raise InvalidInputError(
                    "The field 'end_of_hole' must be a number. Check the input."
                ) from exc
            if not np.isfinite(eoh):
                raise InvalidInputError(
                    "The field 'end_of_hole' must be a finite number. Check the input."
                )
            object.__setattr__(self, "end_of_hole", eoh)
        object.__setattr__(self, "collar", collar)
        object.__setattr__(self, "survey", survey)
        object.__setattr__(self, "intervals", intervals)
        object.__setattr__(
            self, "properties", _properties(self.properties, len(survey), "properties")
        )
        object.__setattr__(
            self,
            "interval_properties",
            _properties(self.interval_properties, len(intervals), "interval_properties"),
        )


SceneObject = Union[Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole]
