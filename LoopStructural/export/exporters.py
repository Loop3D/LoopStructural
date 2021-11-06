"""
Routines to export geological model data to file in a variety of formats
"""
import logging
import os
from pyevtk.hl import unstructuredGridToVTK, pointsToVTK
from pyevtk.vtk import VtkTriangle
import numpy as np
from skimage.measure import marching_cubes

from LoopStructural.utils.helper import create_box
from LoopStructural.export.file_formats import FileFormat


from LoopStructural.utils import getLogger
logger = getLogger(__name__)

def write_feat_surfs(model, file_name, file_format=FileFormat.GOCAD, features=[]):
    """
    Writes out features from a model as 3d surfaces

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path, but without the file extension
    file_format: export.fileformats.FileFormat object
        Desired format of exported file. Supports VTK
    labels: OPTIONAL list of feature names to export, if omitted all faults are exported

    Returns
    -------
    True if successful

    """
    if file_format == FileFormat.GOCAD:
        return _write_feat_surfs_gocad(model, file_name, features)

    logger.warning("Cannot export to surface file - format {} not supported yet".format(str(file_format)))


def _write_feat_surfs_gocad(model, file_name, features=[]):
    for feature in model.features:
        x = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0], model.nsteps[0])
        y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1], model.nsteps[1])
        z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2], model.nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        points = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        val = feature.evaluate_value(points)
        slices_ = [0.0]

        step_vector = np.array([x[1] - x[0], y[1] - y[0], z[1] - z[0]])
        for i, isovalue in enumerate(slices_):
            logger.info("Creating isosurface of {feature.name} at {isovalue}")

            if isovalue > np.nanmax(val) or isovalue < np.nanmin(val):
                logger.warning("For {feature.name} isovalue doesn't exist inside bounding box, skipping")
                continue
            try:
                verts, faces, normals, values = marching_cubes(
                    val.reshape(model.nsteps, order='C'),
                    isovalue,
                    spacing=step_vector)
                verts += np.array([model.bounding_box[0, 0], model.bounding_box[0, 1], model.bounding_box[1, 2]])
                model.rescale(verts)

            except (ValueError, RuntimeError) as e:
                print(e)
                logger.warning(f"Cannot isosurface {feature.name} at {isovalue}, skipping")
                continue
            clean_feature = feature.name.replace(' ','-')
            with open(f"{file_name}_{clean_feature}.TS", "w") as fd:
                fd.write(f"""GOCAD TSurf 1 
HEADER {{
*solid*color: #ffa500
ivolmap: false
imap: false
name: {feature.name}
}}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME Default
PROJECTION Unknown
DATUM Unknown
AXIS_NAME X Y Z
AXIS_UNIT m m m
ZPOSITIVE Elevation
END_ORIGINAL_COORDINATE_SYSTEM
GEOLOGICAL_FEATURE {clean_feature}
GEOLOGICAL_TYPE fault
PROPERTY_CLASS_HEADER X {{
kind: X
unit: m
}}
PROPERTY_CLASS_HEADER Y {{
kind: Y
unit: m
}}
PROPERTY_CLASS_HEADER Z {{
kind: Z
unit: m
is_z: on
}}
PROPERTY_CLASS_HEADER vector3d {{
kind: Length
unit: m
}}
TFACE
""")
                nan_list = []
                v_idx = 1
                v_map = {}
                for idx, vert in enumerate(verts):
                    if not np.isnan(vert[0]) and not np.isnan(vert[1]) and not np.isnan(vert[2]):
                        fd.write(f"VRTX {v_idx:} {vert[0]} {vert[1]} {vert[2]} \n")
                        v_map[idx] = v_idx
                        v_idx += 1
                for face in faces:
                    if face[0] in v_map and face[1] in v_map and face[2] in v_map:
                        fd.write(f"TRGL {v_map[face[0]]} {v_map[face[1]]} {v_map[face[2]]} \n")
                fd.write("END\n")

    return True


def write_cubeface(model, file_name, data_label, nsteps, file_format):
    """
    Writes out the model as a cuboid with six rectangular surfaces

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path, but without the file extension
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions
    file_format: export.fileformats.FileFormat object
        Desired format of exported file. Supports VTK

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
        Name of file that model is exported to, including path, but without the file extension
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions
    file_format: export.fileformats.FileFormat object
        Desired format of exported file. Supports VTK and GOCAD

    Returns
    -------
    True if successful

    """
    if file_format == FileFormat.VTK:
        return _write_vol_evtk(model, file_name, data_label, nsteps)
    if file_format == FileFormat.GOCAD:
        return _write_vol_gocad(model, file_name, data_label, nsteps)

    logger.warning("Cannot export to file - format {} not supported yet".format(str(file_format)))
    return False


def _write_cubeface_evtk(model, file_name, data_label, nsteps, real_coords=True):
    """
    Writes out the model as a cuboid with six rectangular surfaces in VTK unstructured grid format

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path, but without the file extension
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
        unstructuredGridToVTK(file_name, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype,
                              cellData=None, pointData={data_label: val})
    except Exception as e:
        logger.warning("Cannot export cuboid surface to VTK file {}: {}".format(file_name, str(e)))
        return False
    return True


def _write_vol_evtk(model, file_name, data_label, nsteps, real_coords=True):
    """
    Writes out the model as a 3d volume grid in VTK points format

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path, but without the file extension
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
    # xyz is N x 3 vector array
    xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
    vals = model.evaluate_model(xyz, scale=False)
    if real_coords:
        model.rescale(xyz)

    # Define vertices - xyz.shape[0] is length of vector array
    x = np.zeros(xyz.shape[0])
    y = np.zeros(xyz.shape[0])
    z = np.zeros(xyz.shape[0])
    for i in range(xyz.shape[0]):
        x[i], y[i], z[i] = xyz[i][0], xyz[i][1], xyz[i][2]

    # Write to grid
    try:
        pointsToVTK(file_name, x, y, z, data={data_label: vals})
    except Exception as e:
        logger.warning("Cannot export volume to VTK file {}: {}".format(file_name, str(e)))
        return False
    return True


def _write_vol_gocad(model, file_name, data_label, nsteps, real_coords=True):
    """
    Writes out the model as a 3d volume grid in GOCAD VOXET object format

    Parameters
    ----------
    model : GeologicalModel object
        Geological model to export
    file_name : string
        Name of file that model is exported to, including path, but without the file extension
    data_label : string
        A data label to insert into export file
    nsteps : np.array([num-x-steps, num-y-steps, num-z-steps])
        3d array dimensions

    Returns
    -------
    True if successful

    """
    # Define grid spacing in model scale coords
    loop_X = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0], nsteps[0])
    loop_Y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1], nsteps[1])
    loop_Z = np.linspace(model.bounding_box[0, 2], model.bounding_box[1, 2], nsteps[2])

    # Generate model values in 3d grid
    xx, yy, zz = np.meshgrid(loop_X, loop_Y, loop_Z, indexing='ij')
    # xyz is N x 3 vector array
    xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
    vals = model.evaluate_model(xyz, scale=False)
    # Use FORTRAN style indexing for GOCAD VOXET
    vol_vals = np.reshape(vals, nsteps, order='F')
    bbox = model.bounding_box[:]

    # Convert bounding box to real world scale coords
    if real_coords:
        model.rescale(bbox)

    # If integer values
    if type(vals[0]) is np.int64:
        d_type = np.int8
        no_data_val = None
        prop_esize = 1
        prop_storage_type = "Octet"

    # If float values
    elif type(vals[0]) is np.float32:
        d_type = np.dtype('>f4')
        no_data_val = -999999.0
        prop_esize = 4
        prop_storage_type = "Float"
    else:
        logger.warning("Cannot export volume to GOCAD VOXET file: Unsupported type {}".format(type(vals[0])))
        return False

    # Write out VOXET file
    vo_filename = file_name + ".vo"
    data_filename = file_name + "@@"
    try:
        with open(vo_filename, "w") as fp:
            fp.write("""GOCAD Voxet 1
HEADER {{
name: {name}
}}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME Default
AXIS_NAME "X" "Y" "Z"
AXIS_UNIT "m" "m" "m"
ZPOSITIVE Elevation
END_ORIGINAL_COORDINATE_SYSTEM
AXIS_O 0.000000 0.000000 0.000000
AXIS_U 1.000000 0.000000 0.000000
AXIS_V 0.000000 1.000000 0.000000
AXIS_W 0.000000 0.000000 1.000000
AXIS_MIN {axismin1} {axismin2} {axismin3}
AXIS_MAX {axismax1} {axismax2} {axismax3}
AXIS_N {nsteps1} {nsteps2} {nsteps3}
AXIS_NAME "X" "Y" "Z"
AXIS_UNIT "m" "m" "m"
AXIS_TYPE even even even
PROPERTY 1 {propname}
PROPERTY_CLASS 1 {propname}
PROP_UNIT 1 {propname}
PROPERTY_CLASS_HEADER 1 {propname} {{
}}
PROPERTY_SUBCLASS 1 QUANTITY {prop_storage_type}
""".format(name=os.path.basename(file_name),
                nsteps1=nsteps[0], nsteps2=nsteps[1], nsteps3=nsteps[2],
                axismin1=bbox[0, 0], axismin2=bbox[0, 1], axismin3=bbox[0, 2],
                axismax1=bbox[1, 0], axismax2=bbox[1, 1], axismax3=bbox[1, 2],
                propname=data_label, prop_storage_type=prop_storage_type))
            if no_data_val is not None:
                fp.write("PROP_NO_DATA_VALUE 1 {no_data_val}\n".format(no_data_val=no_data_val))
            fp.write("""PROP_ETYPE 1 IEEE
PROP_FORMAT 1 RAW
PROP_ESIZE 1 {prop_esize}
PROP_OFFSET 1 0
PROP_FILE 1 {prop_file}
END\n""".format(prop_file=data_filename, prop_esize=prop_esize))
    except IOError as exc:
        logger.warning("Cannot export volume to GOCAD VOXET file {}: {}".format(vo_filename, str(exc)))
        return False

    # Write out accompanying binary data file
    export_vals = np.array(vol_vals, dtype=d_type)
    try:
        with open(data_filename, "wb") as fp:
            export_vals.tofile(fp)
    except IOError as exc:
        logger.warning("Cannot export volume to GOCAD VOXET data file {}: {}".format(data_filename, str(exc)))
        return False
    return True
