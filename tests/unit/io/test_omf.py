import numpy as np
import pytest

omf = pytest.importorskip("omf")

from LoopStructural.export.omf_wrapper import (
    add_pointset_to_omf,
    add_structured_grid_to_omf,
    add_surface_to_omf,
    get_cell_attributes,
    get_point_attributed,
    get_project,
)
from LoopStructural.geometry import Surface, ValuePoints


class _FakeLoopObject:
    """Minimal duck-typed stand-in for Surface/StructuredGrid used by the
    attribute helpers - only `properties`/`cell_properties` are accessed."""

    def __init__(self, properties=None, cell_properties=None):
        self.properties = properties
        self.cell_properties = cell_properties


def _triangle_surface(name="TestSurface", properties=None):
    return Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        triangles=np.array([[0, 1, 2]]),
        name=name,
        properties=properties,
    )


def test_get_project_returns_new_project_when_file_missing(tmp_path):
    filename = tmp_path / "does_not_exist.omf"
    project = get_project(str(filename))
    assert isinstance(project, omf.Project)
    assert project.name == "LoopStructural Model"
    assert len(project.elements) == 0


def test_get_cell_attributes_empty_when_no_properties():
    obj = _FakeLoopObject(cell_properties=None)
    assert get_cell_attributes(obj) == []


def test_get_cell_attributes_scalar_property():
    obj = _FakeLoopObject(cell_properties={"rock": np.array([1.0, 2.0, 3.0])})
    attributes = get_cell_attributes(obj)
    assert len(attributes) == 1
    assert attributes[0].name == "rock"
    assert attributes[0].location == "faces"
    assert np.allclose(attributes[0].array.array, [1.0, 2.0, 3.0])


def test_get_cell_attributes_multi_column_property_split_by_index():
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    obj = _FakeLoopObject(cell_properties={"vec": values})
    attributes = get_cell_attributes(obj)
    names = sorted(a.name for a in attributes)
    assert names == ["vec_0", "vec_1"]
    for attribute in attributes:
        assert attribute.location == "faces"


def test_get_point_attributed_empty_when_no_properties():
    obj = _FakeLoopObject(properties=None)
    assert get_point_attributed(obj) == []


def test_get_point_attributed_scalar_property():
    obj = _FakeLoopObject(properties={"value": np.array([1.0, 2.0])})
    attributes = get_point_attributed(obj)
    assert len(attributes) == 1
    assert attributes[0].name == "value"
    assert attributes[0].location == "vertices"
    assert np.allclose(attributes[0].array.array, [1.0, 2.0])


def test_add_surface_to_omf_round_trip(tmp_path):
    filename = tmp_path / "surface.omf"
    surf = _triangle_surface(properties={"myprop": np.array([1.0, 2.0, 3.0])})

    add_surface_to_omf(surf, str(filename))
    assert filename.exists()

    project = omf.OMFReader(str(filename)).get_project()
    assert len(project.elements) == 1
    element = project.elements[0]
    assert element.name == "TestSurface"
    assert np.allclose(element.geometry.vertices.array, surf.vertices)
    assert np.array_equal(element.geometry.triangles.array, surf.triangles)
    assert [d.name for d in element.data] == ["myprop"]
    assert np.allclose(element.data[0].array.array, [1.0, 2.0, 3.0])


def test_add_surface_to_omf_appends_to_existing_project_file(tmp_path):
    filename = tmp_path / "two_surfaces.omf"
    add_surface_to_omf(_triangle_surface(name="First"), str(filename))
    add_surface_to_omf(_triangle_surface(name="Second"), str(filename))

    project = omf.OMFReader(str(filename)).get_project()
    assert {element.name for element in project.elements} == {"First", "Second"}


def test_add_structured_grid_to_omf_raises_not_implemented():
    # add_structured_grid_to_omf explicitly rejects structured grids - the
    # real implementation below it is commented out, so structured grids
    # cannot currently be exported to omf.
    with pytest.raises(NotImplementedError, match="cannot store structured grids"):
        add_structured_grid_to_omf(object(), "unused.omf")


@pytest.mark.xfail(
    reason=(
        "add_pointset_to_omf calls omf.PointSetElement(vertices=..., attributes=...) "
        "directly, but the installed omf package (mira-omf) requires "
        "geometry=omf.PointSetGeometry(vertices=...) and data=attributes instead. "
        "This raises AttributeError: 'Keyword input is not a known property of "
        "PointSetElement' - a bug in LoopStructural/export/omf_wrapper.py."
    ),
    strict=True,
    raises=AttributeError,
)
def test_add_pointset_to_omf_round_trip(tmp_path):
    filename = tmp_path / "points.omf"
    points = ValuePoints(
        locations=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
        values=np.array([10.0, 20.0, 30.0]),
        name="TestPoints",
    )

    add_pointset_to_omf(points, str(filename))

    project = omf.OMFReader(str(filename)).get_project()
    assert project.elements[0].name == "TestPoints"
