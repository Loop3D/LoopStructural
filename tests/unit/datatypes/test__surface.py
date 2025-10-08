import numpy as np
import pytest
from LoopStructural.datatypes._surface import Surface


def test_surface_creation():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(
        vertices=vertices, triangles=triangles, normals=normals, name=name, values=values
    )

    assert np.array_equal(surface.vertices, vertices)
    assert np.array_equal(surface.triangles, triangles)
    assert np.array_equal(surface.normals, normals)
    assert surface.name == name
    assert np.array_equal(surface.values, values)


def test_surface_vtk():
    import importlib.util

    if importlib.util.find_spec('pyvista') is None:
        pytest.skip("pyvista is required for vtk support")
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(
        vertices=vertices, triangles=triangles, normals=normals, name=name, values=values
    )
    vtk_surface = surface.vtk()

    assert vtk_surface.n_points == len(vertices)
    assert vtk_surface.n_cells == len(triangles)
    assert np.array_equal(vtk_surface.point_data["values"], values)


def test_surface_to_dict():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(
        vertices=vertices, triangles=triangles, normals=normals, name=name, values=values
    )
    surface_dict = surface.to_dict()

    assert np.array_equal(surface_dict["vertices"], vertices)
    assert np.array_equal(surface_dict["triangles"], triangles)
    assert np.array_equal(surface_dict["normals"], normals)
    assert surface_dict["name"] == name
    assert np.array_equal(surface_dict["values"], values)

def test_surface_remove_nan_vertices():
    # Vertices 1 and 3 are NaN, should be removed
    vertices = np.array([
        [0, 0, 0],
        [np.nan, np.nan, np.nan],
        [1, 0, 0],
        [np.nan, 1, 0],
        [0, 1, 0]
    ])
    triangles = np.array([
        [0, 2, 4],  # valid
        [1, 2, 4],  # contains NaN vertex 1, should be removed
        [0, 3, 4]   # contains NaN vertex 3, should be removed
    ])
    normals = np.array([
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1]
    ])
    values = np.array([10, 20, 30, 40, 50])
    properties = {"prop": np.array([1, 2, 3, 4, 5])}
    cell_properties = {"cell_prop": np.array([100, 200, 300])}

    surface = Surface(
        vertices=vertices,
        triangles=triangles,
        normals=normals,
        values=values,
        properties=properties,
        cell_properties=cell_properties
    )

    # Only valid vertices and triangles should remain
    assert surface.vertices.shape[0] == 3
    assert surface.triangles.shape[0] == 1
    # The only triangle left should reference the new indices [0,1,2]
    assert np.array_equal(surface.triangles[0], [0, 1, 2])
    # Normals, values, and properties should be filtered accordingly
    assert surface.normals.shape[0] == 3
    assert np.array_equal(surface.values, [10, 30, 50])
    assert np.array_equal(surface.properties["prop"], [1, 3, 5])
    assert np.array_equal(surface.cell_properties["cell_prop"], [100])
