"""
Routines to export geological model data to file in a variety of formats
"""
import logging
from pyevtk.hl import unstructuredGridToVTK, pointsToVTK
from pyevtk.vtk import VtkTriangle
import numpy as np

from LoopStructural.utils.helper import create_box
from LoopStructural.export.file_formats import FileFormat

 
logger = logging.getLogger(__name__)


def write_cubeface(model, file_name, data_label, nsteps, file_format):
    """
    Writes out the model as a cuboid with six rectangular surfaces

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions
    file_format: export.fileformats.FileFormat object
        Desired format of exported file

    Returns
    -------
    True if successful

    """
    if file_format == FileFormat.VTK:
        return _write_cubeface_evtk(model, file_name, data_label, nsteps)

    logger.warning("Cannot export to file - format {} not supported yet".format(str(file_format)))
    return False


def write_vol(model, file_name, data_label, nsteps, file_format):
    """
    Writes out the model as a 3d volume grid

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions
    file_format: export.fileformats.FileFormat object
        Desired format of exported file

    Returns
    -------
    True if successful

    """
    if file_format == FileFormat.VTK:
        return _write_vol_evtk(model, file_name, data_label, nsteps)

    logger.warning("Cannot export to file - format {} not supported yet".format(str(file_format)))
    return False


def _write_cubeface_evtk(model, file_name, data_label, nsteps, real_coords=True):
    """
    Writes out the model as a cuboid with six rectangular surfaces

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions

    Returns
    -------
    True if successful

    """
    # Evaluate model at points
    points, tri = create_box(model.bounding_box, nsteps)
    val = model.evaluate_model(points, scale=False)
    if real_coords:
        model.rescale(points)

    # Define vertices
    x = np.zeros(points.shape[0])
    y = np.zeros(points.shape[0])
    z = np.zeros(points.shape[0])
    for i in range(points.shape[0]):
        x[i], y[i], z[i] = points[i][0], points[i][1], points[i][2]

    # Define connectivity or vertices that belongs to each element
    conn = np.zeros(tri.shape[0] * 3)
    for i in range(tri.shape[0]):
        conn[i*3], conn[i*3+1], conn[i*3+2] = tri[i][0], tri[i][1], tri[i][2]

    # Define offset of last vertex of each element
    offset = np.zeros(tri.shape[0])
    for i in range(tri.shape[0]):
        offset[i] = (i+1)*3

    # Define cell types
    ctype = np.full(tri.shape[0], VtkTriangle.tid)

    try:
        unstructuredGridToVTK(file_name, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, cellData = None, pointData = {data_label: val})
    except Exception as e:
        logger.warning("Cannot export cuboid surface to file {}: {}".format(file_name, str(e)))
        return False
    return True 


def _write_vol_evtk(model, file_name, data_label, nsteps, real_coords=True):
    """
    Writes out the model as a 3d volume grid

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions

    Returns
    -------
    True if successful

    """
    # Define grid spacing
    loop_X = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0], nsteps[0])
    loop_Y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1], nsteps[1])
    loop_Z = np.linspace(model.bounding_box[0, 2], model.bounding_box[1, 2], nsteps[2])

    # Generate model values in 3d grid
    xx, yy, zz = np.meshgrid(loop_X, loop_Y, loop_Z, indexing='ij')
    xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
    vals = model.evaluate_model(xyz, scale=False)
    if real_coords:
        model.rescale(xyz)

    # Define vertices
    x = np.zeros(xyz.shape[0])
    y = np.zeros(xyz.shape[0])
    z = np.zeros(xyz.shape[0])
    for i in range(xyz.shape[0]):
        x[i], y[i], z[i] = xyz[i][0], xyz[i][1], xyz[i][2]

    # Write to grid
    try:
        pointsToVTK(file_name, x, y, z, data= { data_label: vals})
    except Exception as e:
        logger.warning("Cannot export volume to file {}: {}".format(file_name, str(e)))
        return False
    return True 


