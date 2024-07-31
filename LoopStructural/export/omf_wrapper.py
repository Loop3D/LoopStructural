try:
    import omf
except ImportError:
    raise ImportError(
        "You need to install the omf package to use this feature. "
        "You can install it with: pip install --pre omf"
    )
import numpy as np


def get_project(filename):
    try:
        project = omf.load(filename)
    except FileNotFoundError:
        project = omf.Project(name='LoopStructural Model')
    return project


def get_cell_attributes(loopobject):
    attributes = []
    if loopobject.cell_properties:
        attributes += [
            omf.NumericAttribute(name=k, array=v, location="faces")
            for k, v in loopobject.cell_properties.items()
        ]
    return attributes


def get_point_attributed(loopobject):
    attributes = []
    if loopobject.properties:
        attributes += [
            omf.NumericAttribute(name=k, array=v, location="vertices")
            for k, v in loopobject.properties.items()
        ]
    return attributes


def add_surface_to_omf(surface, filename):

    attributes = []
    attributes += get_cell_attributes(surface)
    attributes += get_point_attributed(surface)
    surface = omf.Surface(
        vertices=surface.vertices,
        triangles=surface.triangles,
        attributes=attributes,
        name=surface.name,
    )
    project = get_project(filename)

    project.elements += [surface]
    omf.save(project, filename, mode='w')


def add_pointset_to_omf(points, filename):

    attributes = []
    attributes += get_point_attributed(points)

    points = omf.PointSet(
        vertices=points.locations,
        attributes=attributes,
        name=points.name,
    )

    project = get_project(filename)
    project.elements += [points]
    omf.save(project, filename, mode='w')


def add_structured_grid_to_omf(grid, filename):
    print('Open Mining Format cannot store structured grids')
    return
    # attributes = []
    # attributes += get_cell_attributes(grid)
    # attributes += get_point_attributed(grid)

    # vol = omf.TensorGridBlockModel(
    #     name=grid.name,
    #     tensor_u=np.ones(grid.nsteps[0]).astype(float),
    #     tensor_v=np.ones(grid.nsteps[0]).astype(float),
    #     tensor_w=np.ones(grid.nsteps[0]).astype(float),
    #     origin=grid.origin,
    #     attributes=attributes,
    # )
    # project = get_project(filename)
    # project.elements += [vol]
    # omf.save(project, filename, mode='w')
