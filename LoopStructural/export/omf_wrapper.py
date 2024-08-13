try:
    import omf
except ImportError:
    raise ImportError(
        "You need to install the omf package to use this feature. "
        "You can install it with: pip install mira-omf"
    )
import numpy as np
import datetime
import os


def get_project(filename):
    if os.path.exists(filename):

        try:
            reader = omf.OMFReader(filename)
            project = reader.get_project()
        except (FileNotFoundError, ValueError):
            project = omf.Project(name='LoopStructural Model')
        return project
    else:
        return omf.Project(name='LoopStructural Model')


def get_cell_attributes(loopobject):
    attributes = []
    if loopobject.cell_properties:
        for k, v in loopobject.cell_properties.items():
            v = np.array(v)
            if len(v.shape) > 1 and v.shape[1] > 1:
                for i in range(v.shape[1]):
                    attributes.append(
                        omf.ScalarData(name=f'{k}_{i}', array=v[:, i], location="faces")
                    )
            else:
                attributes.append(omf.ScalarData(name=k, array=v, location="faces"))

    return attributes


def get_point_attributed(loopobject):
    attributes = []
    if loopobject.properties:
        for k, v in loopobject.properties.items():
            v = np.array(v)
            if len(v.shape) > 1 and v.shape[1] > 1:
                for i in range(v.shape[1]):
                    attributes.append(
                        omf.ScalarData(name=f'{k}_{i}', array=v[:, i], location="vertices")
                    )
            else:
                attributes.append(omf.ScalarData(name=k, array=v, location="vertices"))

    return attributes


def add_surface_to_omf(surface, filename):

    attributes = []
    attributes += get_cell_attributes(surface)
    attributes += get_point_attributed(surface)
    surface = omf.SurfaceElement(
        geometry=omf.SurfaceGeometry(
            vertices=surface.vertices,
            triangles=surface.triangles,
        ),
        data=attributes,
        name=surface.name,
    )
    project = get_project(filename)

    project.elements += [surface]
    project.metadata = {
        "coordinate_reference_system": "epsg 3857",
        "date_created": datetime.datetime.utcnow(),
        "version": "v1.3",
        "revision": "10",
    }
    omf.OMFWriter(project, filename)


def add_pointset_to_omf(points, filename):

    attributes = []
    attributes += get_point_attributed(points)

    points = omf.PointSetElement(
        vertices=points.locations,
        attributes=attributes,
        name=points.name,
    )

    project = get_project(filename)
    project.elements += [points]
    omf.OMFWriter(project, filename)


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
