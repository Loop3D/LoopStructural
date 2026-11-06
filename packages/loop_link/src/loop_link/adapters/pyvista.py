"""Adapter that turns pyvista meshes and grids into loop_link scene objects."""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np

from loop_link.errors import InvalidInputError, MissingExtraError
from loop_link.objects import (
    BlockModel,
    Lines,
    Points,
    SceneObject,
    Surface,
    TetraMesh,
    Vectors,
)

try:
    import pyvista as _pv
except ImportError:  # pragma: no cover - depends on the environment
    _pv = None

_DRILLHOLE_TEXT = (
    "A pyvista mesh cannot hold a drillhole. Use loop_link.Drillhole or "
    "loop_link.adapters.pandas.drillholes_from_pandas."
)

# Arrays that pyvista adds by itself. The adapter skips them.
_SKIPPED_ARRAYS = ("Normals", "Texture Coordinates", "TextureCoordinates")


def _refuse(reason: str) -> InvalidInputError:
    """Make an error that gives the reason and the accepted inputs."""
    from loop_link.adapters import _ACCEPTED_INPUTS

    return InvalidInputError(f"{reason} {_ACCEPTED_INPUTS}")


def _collect(data: Any, where: str) -> dict[str, np.ndarray]:  # noqa: ANN401
    """Copy the arrays of a pyvista data set to a dictionary of numpy arrays.

    Parameters
    ----------
    data : Any
        A pyvista ``point_data`` or ``cell_data`` object.
    where : str
        The name of the data set, used in error text.

    Returns
    -------
    dict of str to numpy.ndarray
        The 1-D arrays, without the arrays that pyvista adds itself.

    Raises
    ------
    InvalidInputError
        An array is not 1-D.
    """
    result: dict[str, np.ndarray] = {}
    for name in list(data.keys()):
        if name in _SKIPPED_ARRAYS:
            continue
        array = np.asarray(data[name])
        if array.ndim != 1:
            raise InvalidInputError(
                f"The {where} array '{name}' has {array.ndim} dimensions. "
                "Only 1-D arrays are accepted. Split the array into 1-D arrays."
            )
        result[name] = array
    return result


def _polydata_lines(mesh: Any) -> tuple[np.ndarray, np.ndarray]:  # noqa: ANN401
    """Return the segments of a PolyData and the source line of each segment."""
    flat = np.asarray(mesh.lines)
    segments: list[list[int]] = []
    owners: list[int] = []
    position = 0
    line = 0
    while position < len(flat):
        count = int(flat[position])
        ids = flat[position + 1 : position + 1 + count]
        for first, second in itertools.pairwise(ids):
            segments.append([int(first), int(second)])
            owners.append(line)
        position += count + 1
        line += 1
    return np.array(segments, dtype=np.int64).reshape(-1, 2), np.array(owners, dtype=np.int64)


def _convert_polydata(
    mesh: Any,  # noqa: ANN401
    vectors: str | None,
    vectors_kind: str,
) -> SceneObject:
    """Convert a PolyData to Lines, Surface, Points or Vectors."""
    if mesh.n_lines > 0:
        segments, owners = _polydata_lines(mesh)
        # Cell data is ordered: vertices, lines, polygons, strips.
        first = int(mesh.n_verts)
        cell_properties = {
            name: values[first : first + int(mesh.n_lines)][owners]
            for name, values in _collect(mesh.cell_data, "cell").items()
        }
        return Lines(
            vertices=np.asarray(mesh.points),
            segments=segments,
            properties=_collect(mesh.point_data, "point"),
            cell_properties=cell_properties,
        )
    if len(mesh.faces) > 0:
        if mesh.n_verts > 0 or mesh.n_strips > 0:
            raise InvalidInputError(
                "The PolyData mixes polygons with vertices or strips. "
                "Split the mesh into one mesh for each kind of cell."
            )
        if not mesh.is_all_triangles:
            mesh = mesh.triangulate()
        return Surface(
            vertices=np.asarray(mesh.points),
            triangles=np.asarray(mesh.faces).reshape(-1, 4)[:, 1:],
            properties=_collect(mesh.point_data, "point"),
            cell_properties=_collect(mesh.cell_data, "cell"),
        )
    points = np.asarray(mesh.points)
    if vectors is None:
        return Points(vertices=points, properties=_collect(mesh.point_data, "point"))
    if vectors not in mesh.point_data:
        raise InvalidInputError(
            f"The point array '{vectors}' does not exist. Give the name of a point array."
        )
    directions = np.asarray(mesh.point_data[vectors])
    properties = _collect(_Without(mesh.point_data, vectors), "point")
    return Vectors.from_directions(points, directions, kind=vectors_kind, properties=properties)


class _Without:
    """A read-only view of a pyvista data set that hides one array."""

    def __init__(self, data: Any, hidden: str) -> None:  # noqa: ANN401
        self._data = data
        self._hidden = hidden

    def keys(self) -> list[str]:
        """Return the array names without the hidden name."""
        return [name for name in self._data if name != self._hidden]

    def __getitem__(self, name: str) -> Any:  # noqa: ANN401
        return self._data[name]


def _check_direction(mesh: Any) -> None:  # noqa: ANN401
    """Raise an error when an ImageData has a rotation."""
    direction = getattr(mesh, "direction_matrix", None)
    if direction is not None and not np.allclose(np.asarray(direction), np.eye(3)):
        raise InvalidInputError(
            "The image data has a rotation that the viewer cannot show. "
            "Resample the image data to a grid without rotation, or remove the rotation "
            "(set direction_matrix to the identity matrix), then try again."
        )


def _even_step(axis: np.ndarray, label: str) -> float:
    """Return the step of an axis, or raise an error when the steps are not equal."""
    if len(axis) < 2:
        raise _refuse(f"The {label} axis of the rectilinear grid has less than 2 nodes.")
    step = float(axis[1] - axis[0])
    if not np.allclose(np.diff(axis), step):
        raise _refuse(
            f"The {label} axis of the rectilinear grid has unequal steps. "
            "A BlockModel needs equal steps on each axis."
        )
    return step


def _convert_block(mesh: Any, pv: Any) -> BlockModel:  # noqa: ANN401
    """Convert an ImageData or a RectilinearGrid to a BlockModel."""
    if isinstance(mesh, pv.ImageData):
        _check_direction(mesh)
        origin = np.asarray(mesh.origin, dtype=np.float64)
        sizes = np.asarray(mesh.spacing, dtype=np.float64)
        cells = np.array(mesh.dimensions) - 1
    else:
        axes = [np.asarray(mesh.x), np.asarray(mesh.y), np.asarray(mesh.z)]
        sizes = np.array([_even_step(a, n) for a, n in zip(axes, "xyz")], dtype=np.float64)
        origin = np.array([a[0] for a in axes], dtype=np.float64)
        cells = np.array([len(a) - 1 for a in axes])
    return BlockModel(
        origin=origin,
        cell_sizes=sizes,
        cells=cells,
        properties=_collect(mesh.point_data, "point"),
        cell_properties=_collect(mesh.cell_data, "cell"),
    )


def from_pyvista(
    mesh: Any,  # noqa: ANN401
    *,
    kind: str | None = None,
    vectors: str | None = None,
    vectors_kind: str = "normal",
) -> SceneObject:
    """Convert a pyvista mesh or grid to a loop_link scene object.

    Parameters
    ----------
    mesh : pyvista.DataSet
        A ``PolyData``, ``ImageData``, ``RectilinearGrid`` or ``UnstructuredGrid``.
    kind : str, optional
        The name of the scene object type that the caller expects. When it does not
        match the type that the mesh gives, an error is raised.
    vectors : str, optional
        The name of a point array that holds direction vectors. It is used for a
        ``PolyData`` of vertices only, and gives ``Vectors``.
    vectors_kind : str, default "normal"
        How to read the vectors: ``"normal"`` or ``"lineation"``.

    Returns
    -------
    SceneObject
        The converted object. Point data gives ``properties``. Cell data gives
        ``cell_properties``. The order of image data is x fastest, and is not changed.

    Raises
    ------
    MissingExtraError
        pyvista is not installed.
    InvalidInputError
        The mesh type is not accepted, ``kind`` does not match, or an array is not 1-D.
    """
    if _pv is None:
        raise MissingExtraError("pyvista")
    pv = _pv
    if kind == "Drillhole":
        raise InvalidInputError(_DRILLHOLE_TEXT)
    result: SceneObject
    if isinstance(mesh, pv.PolyData):
        result = _convert_polydata(mesh, vectors, vectors_kind)
    elif isinstance(mesh, (pv.ImageData, pv.RectilinearGrid)):
        result = _convert_block(mesh, pv)
    elif isinstance(mesh, pv.UnstructuredGrid):
        types = np.asarray(mesh.celltypes)
        if len(types) == 0 or not np.all(types == pv.CellType.TETRA):
            raise _refuse("An UnstructuredGrid is accepted only when all cells are tetrahedra.")
        result = TetraMesh(
            vertices=np.asarray(mesh.points),
            tetrahedra=np.asarray(mesh.cells_dict[pv.CellType.TETRA]),
            properties=_collect(mesh.point_data, "point"),
            cell_properties=_collect(mesh.cell_data, "cell"),
        )
    else:
        raise _refuse(f"The pyvista type '{type(mesh).__name__}' is not accepted.")
    if kind is not None and kind != type(result).__name__:
        raise _refuse(
            f"The kind '{kind}' does not match the mesh. The mesh gives '{type(result).__name__}'."
        )
    return result
