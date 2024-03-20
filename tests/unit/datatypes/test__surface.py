import numpy as np
import pytest
from LoopStructural.datatypes._surface import Surface


def test_surface_creation():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(vertices, triangles, normals, name, values)

    assert np.array_equal(surface.vertices, vertices)
    assert np.array_equal(surface.triangles, triangles)
    assert np.array_equal(surface.normals, normals)
    assert surface.name == name
    assert np.array_equal(surface.values, values)


def test_surface_vtk():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(vertices, triangles, normals, name, values)
    vtk_surface = surface.vtk

    assert vtk_surface.n_points == len(vertices)
    assert vtk_surface.n_cells == len(triangles)
    assert np.array_equal(vtk_surface.point_data["values"], values)


def test_surface_to_dict():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])
    normals = np.array([[0, 0, 1]])
    name = "surface1"
    values = np.array([1, 2, 3])

    surface = Surface(vertices, triangles, normals, name, values)
    surface_dict = surface.to_dict()

    assert np.array_equal(surface_dict["vertices"], vertices)
    assert np.array_equal(surface_dict["triangles"], triangles)
    assert np.array_equal(surface_dict["normals"], normals)
    assert surface_dict["name"] == name
    assert np.array_equal(surface_dict["values"], values)
