"""
Routines to export geological model data to file in a variety of formats
"""

import os
from pyevtk.hl import unstructuredGridToVTK, pointsToVTK
from pyevtk.vtk import VtkTriangle
import numpy as np
from skimage.measure import marching_cubes

from LoopStructural.utils.helper import create_box
from LoopStructural.export.file_formats import FileFormat
from LoopStructural.datatypes import Surface

from ..utils import getLogger

logger = getLogger(__name__)


def write_feat_surfs(
    model, featurename, file_format=FileFormat.NUMPY, file_name=None, isovalue=0.0
):
    """
    Writes out features from a model as 3d surfaces

    Parameters
    ----------
    model: GeologicalModel object
        Geological model to export
    file_name: string
        Name of file that model is exported to, including path, but without the file extension
    file_format: export.fileformats.FileFormat object
        OPTIONAL desired format of exported file. Supports GOCAD, VTK & NUMPY. Default is NUMPY
    target_feats: list
        OPTIONAL list of feature names to export, if omitted all faults are exported
    isovalue: float
        OPTIONAL isovalue point at which surface is generated, default is 0.0

    Returns
    -------
    Tuple of (boolean, [ SimpleNamespace() ...  ])
    If successful, boolean is True and SimpleNamespace() objects have the following attributes:
        verts: vertices, numpy ndarray with dtype = float64 & shape = (N,3)
        faces: faces, numpy ndarray with dtype = int32 & shape = (M,3)
        values: values, numpy ndarray with dtype = float32 & shape = (N,)
        normals: normals, numpy ndarray with dtype = float32 & shape = (N,3)
        name: name of feature e.g. fault or supergroup, string
    If not successful, it returns (False, [])
    """
    # Set up the type of file writer function
    if file_format == FileFormat.GOCAD:
        write_fn = _write_feat_surfs_gocad
    elif file_format == FileFormat.VTK:
        write_fn = _write_feat_surfs_evtk
    elif file_format == FileFormat.NUMPY:
        write_fn = None
    else:
        logger.warning(f"Cannot export to surface file - format {file_format} not supported yet")
        return False, []
    # Skip if not a requested feature
    if featurename not in model:
        logger.warning("{featurename} is not in the model, skipping")
        return False, []
    has_feats = True

    x = np.linspace(model.bounding_box.bb[0, 0], model.bounding_box.bb[1, 0], model.nsteps[0])
    y = np.linspace(model.bounding_box.bb[0, 1], model.bounding_box.bb[1, 1], model.nsteps[1])
    z = np.linspace(model.bounding_box.bb[1, 2], model.bounding_box.bb[0, 2], model.nsteps[2])
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    points = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
    val = model[featurename].evaluate_value(points)
    step_vector = np.array([x[1] - x[0], y[1] - y[0], z[1] - z[0]])
    logger.info(f"Creating isosurface of {featurename} at {isovalue}")

    if isovalue > np.nanmax(val) or isovalue < np.nanmin(val):
        logger.warning(
            f"For {featurename} isovalue {isovalue} doesn't exist inside bounding box, skipping"
        )
        return False, []
    try:
        verts, faces, normals, values = marching_cubes(
            val.reshape(model.nsteps, order="C"), isovalue, spacing=step_vector
        )
        verts += np.array(
            [
                model.bounding_box.bb[0, 0],
                model.bounding_box.bb[0, 1],
                model.bounding_box.bb[1, 2],
            ]
        )
        model.rescale(verts)
        surf = Surface(
            vertices=verts,
            triangles=faces,
            normals=normals,
            values=values,
            name=featurename.replace(" ", "-"),
        )

    except (ValueError, RuntimeError) as e:
        logger.debug(f"Exception creating feature surface {featurename}: {e}")
        logger.warning(f"Cannot isosurface {featurename} at {isovalue}, skipping")
        return False, []

    if not has_feats:
        logger.warning("Cannot locate features in model")
        return False, []

    # Call the file writer function
    result = True
    if write_fn is not None:
        result = write_fn(surf, file_name)
    return result


def _write_feat_surfs_evtk(surf, file_name):
    """
    Writes out an unstructured VTK file containing a 2d fault surface

    Parameters
    ----------
    surf_list: [ SimpleNamespace() ... ]
        Details of the surfaces, as a list of SimpleNamespace() objects. Fields are:
            verts: vertices, numpy ndarray with dtype = float64 & shape = (N,3)
            faces: faces, numpy ndarray with dtype = int32 & shape = (M,3)
            values: values, numpy ndarray with dtype = float32 & shape = (N,)
            normals: normals, numpy ndarray with dtype = float32 & shape = (N,3)
            name: name of feature e.g. fault or supergroup, string

    file_name: string
        file name to be written out without extension

    Returns
    -------
    True if successful

    """
    connectivity = np.zeros(0)
    offsets = np.zeros(0)
    cell_types = np.zeros(0)
    x = np.zeros(0)
    y = np.zeros(0)
    z = np.zeros(0)
    conn_idx = 0
    offset_idx = 0
    pointData = np.zeros(0)

    # Concatenate the geometry of all features into x, y, z, connectivity, offsets & cell_types

    # Accumulate x,y,z values in 1-d array
    x = np.append(x, np.ascontiguousarray(surf.verts[:, 0]))
    y = np.append(y, np.ascontiguousarray(surf.verts[:, 1]))
    z = np.append(z, np.ascontiguousarray(surf.verts[:, 2]))

    # Convert faces to int64
    conn = np.array(surf.faces, dtype=np.int64)

    # Reshape connections to 1d
    conn = conn.reshape(conn.size)
    conn += conn_idx

    # Make offsets into connection array
    offs = np.array(list(range(3 + offset_idx, len(conn) + offset_idx, 3)), dtype=np.int64)

    # Set VTK datatype as triangles
    ctype = np.zeros(offs.size)
    ctype.fill(VtkTriangle.tid)

    # Accumulate values in arrays
    connectivity = np.append(connectivity, conn)

    offsets = np.append(offsets, offs)
    cell_types = np.append(cell_types, ctype)

    pointData = np.append(pointData, surf.values.reshape(surf.values.size))

    # Enable connections to point to next set of vertices
    conn_idx = x.size

    # Enable offsets to point to next set of connections
    offset_idx = connectivity.size

    # Write out file
    try:
        logger.info(f"Writing file {file_name}.vtu")
        unstructuredGridToVTK(
            f"{file_name}",
            x,
            y,
            z,
            connectivity=connectivity,
            offsets=offsets,
            cell_types=cell_types,
            pointData={"values": pointData},
        )
    except Exception as e:
        logger.warning(f"Cannot export fault surface to VTK file {file_name}: {e}")
        return False

    return True


def _write_feat_surfs_gocad(surf, file_name):
    """
    Writes out a GOCAD TSURF file for each surface in list

    Parameters
    ----------
    surf_list: [ SimpleNamespace() ... ]
        Details of the surfaces, as a list of SimpleNamespace() objects. Fields are:
            verts: vertices, numpy ndarray with dtype = float64 & shape = (N,3)
            faces: faces, numpy ndarray with dtype = int32 & shape = (M,3)
            values: values, numpy ndarray with dtype = float32 & shape = (N,)
            normals: normals, numpy ndarray with dtype = float32 & shape = (N,3)
            name: name of feature e.g. fault or supergroup, string

    file_name: string
        Desired filename

    Returns
    -------
    True if successful

    """
    from pathlib import Path

    file_name = Path(file_name).with_suffix(".ts")
    with open(f"{file_name}", "w") as fd:
        fd.write(
            f"""GOCAD TSurf 1 
HEADER {{
*solid*color: #ffa500
ivolmap: false
imap: false
name: {surf.name}
}}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME Default
PROJECTION Unknown
DATUM Unknown
AXIS_NAME X Y Z
AXIS_UNIT m m m
ZPOSITIVE Elevation
END_ORIGINAL_COORDINATE_SYSTEM
GEOLOGICAL_FEATURE {surf.name}
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
"""
        )
        v_idx = 1
        v_map = {}
        for idx, vert in enumerate(surf.vertices):
            if not np.isnan(vert[0]) and not np.isnan(vert[1]) and not np.isnan(vert[2]):
                fd.write(f"VRTX {v_idx:} {vert[0]} {vert[1]} {vert[2]} \n")
                v_map[idx] = v_idx
                v_idx += 1
        for face in surf.triangles:
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

    logger.warning(f"Cannot export to file - format {file_format} not supported yet")
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

    logger.warning(f"Cannot export to file - format {file_format} not supported yet")
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
        conn[i * 3], conn[i * 3 + 1], conn[i * 3 + 2] = tri[i][0], tri[i][1], tri[i][2]

    # Define offset of last vertex of each element
    offset = np.zeros(tri.shape[0])
    for i in range(tri.shape[0]):
        offset[i] = (i + 1) * 3

    # Define cell types
    ctype = np.full(tri.shape[0], VtkTriangle.tid)

    try:
        unstructuredGridToVTK(
            file_name,
            x,
            y,
            z,
            connectivity=conn,
            offsets=offset,
            cell_types=ctype,
            cellData=None,
            pointData={data_label: val},
        )
    except Exception as e:
        logger.warning(f"Cannot export cuboid surface to VTK file {file_name}: {e}")
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
    xyz = model.bounding_box.regular_grid(nsteps=nsteps)
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
        logger.warning(f"Cannot export volume to VTK file {file_name}: {e}")
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
    xyz = model.bounding_box.regular_grid(nsteps=nsteps)

    vals = model.evaluate_model(xyz, scale=False)
    # Use FORTRAN style indexing for GOCAD VOXET
    vol_vals = np.reshape(vals, nsteps, order="F")
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
        d_type = np.dtype(">f4")
        no_data_val = -999999.0
        prop_esize = 4
        prop_storage_type = "Float"
    else:
        logger.warning(
            f"Cannot export volume to GOCAD VOXET file: Unsupported type {type(vals[0])}"
        )
        return False

    # Write out VOXET file
    vo_filename = file_name + ".vo"
    data_filename = file_name + "@@"
    try:
        with open(vo_filename, "w") as fp:
            fp.write(
                f"""GOCAD Voxet 1
HEADER {{
name: {os.path.basename(file_name)}
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
AXIS_MIN {bbox[0, 0]} {bbox[0, 1]} {bbox[0, 2]}
AXIS_MAX {bbox[1, 0]} {bbox[1, 1]} {bbox[1, 2]}
AXIS_N {nsteps[0]} {nsteps[1]} {nsteps[2]}
AXIS_NAME "X" "Y" "Z"
AXIS_UNIT "m" "m" "m"
AXIS_TYPE even even even
PROPERTY 1 {data_label}
PROPERTY_CLASS 1 {data_label}
PROP_UNIT 1 {data_label}
PROPERTY_CLASS_HEADER 1 {data_label} {{
}}
PROPERTY_SUBCLASS 1 QUANTITY {prop_storage_type}
"""
            )
            if no_data_val is not None:
                fp.write(f"PROP_NO_DATA_VALUE 1 {no_data_val}\n")
            fp.write(
                f"""PROP_ETYPE 1 IEEE
PROP_FORMAT 1 RAW
PROP_ESIZE 1 {prop_esize}
PROP_OFFSET 1 0
PROP_FILE 1 {data_filename}
END\n"""
            )
    except IOError as exc:
        logger.warning(f"Cannot export volume to GOCAD VOXET file {vo_filename}: {exc}")
        return False

    # Write out accompanying binary data file
    export_vals = np.array(vol_vals, dtype=d_type)
    try:
        with open(data_filename, "wb") as fp:
            export_vals.tofile(fp)
    except IOError as exc:
        logger.warning(f"Cannot export volume to GOCAD VOXET data file {data_filename}: {exc}")
        return False
    return True
