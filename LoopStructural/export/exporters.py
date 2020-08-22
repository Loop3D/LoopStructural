import logging
import meshio
import numpy as np

from LoopStructural.utils.helper import create_box
from LoopStructural.export.file_formats import FileFormat

 
logger = logging.getLogger(__name__)


def write_cubeface_vol(model, file_name, nsteps, file_format):
    """ Writes out a GeologicalModel as a 6-faced cuboid volume to a file

    Parameters
    ----------
    model : model object
       'GeologicalModel' object
    file_name : name of file
       Name of file that model is exported to, including path
    nsteps : number of steps in volume
       NumPy 3x1 array specifying the point density of the cube in three dimensions
    file_format : exported file format
       An 'export.file_formats.FileFormat' object e.g. 'FileFormat.OBJ'

    Returns
    -------
    True if successful
    """
    points, tri = create_box(model.bounding_box, nsteps)
    model.rescale(points)

    val = model.evaluate_model(points, scale=True)

    # Use meshio to output cube faces
    if file_format in [FileFormat.OBJ, FileFormat.VTK]:

        # Remove 'FileFormat.'
        out_fmt = str(file_format)[11:].lower()

        # Write out mesh
        cells = [("triangle", tri)] 
        try:
            meshio.write_points_cells(file_name, points, cells, file_format=out_fmt)
        except Exception as e:
            logger.error("Cannot export file {}: {}".format(file_name, str(e)))
            return False 

    # Use pyassimp to output cube faces
    elif file_format == FileFormat.GLTF:
        # Not supported yet
        logger.warning("Cannot export to file - GLTF not supported yet")
        return False

    # Unknown file type
    else:
        logger.error("Cannot export to file - unknown FileType")
        return False
    return True
