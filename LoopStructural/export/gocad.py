import re
from pathlib import Path

import numpy as np

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


def _normalise_voxet_property(values, property_name, nsteps):
    array = np.asarray(values)
    expected_shape = tuple(int(step) for step in nsteps)
    expected_size = int(np.prod(expected_shape))

    if array.shape == expected_shape:
        flat_values = array.reshape(-1, order="F")
    else:
        flat_values = np.squeeze(array)
        if flat_values.shape == expected_shape:
            flat_values = flat_values.reshape(-1, order="F")
        elif flat_values.ndim == 1 and flat_values.size == expected_size:
            flat_values = flat_values
        else:
            raise ValueError(
                f"Property '{property_name}' must have shape {expected_shape} or size {expected_size}"
            )

    if np.issubdtype(flat_values.dtype, np.integer):
        if flat_values.size == 0:
            export_dtype = np.int8
            storage_type = "Octet"
            element_size = 1
        elif flat_values.min() >= np.iinfo(np.int8).min and flat_values.max() <= np.iinfo(np.int8).max:
            export_dtype = np.int8
            storage_type = "Octet"
            element_size = 1
        else:
            export_dtype = np.dtype(">i4")
            storage_type = "Integer"
            element_size = 4
        no_data_value = None
    elif np.issubdtype(flat_values.dtype, np.floating):
        export_dtype = np.dtype(">f4")
        storage_type = "Float"
        element_size = 4
        no_data_value = -999999.0
        flat_values = np.nan_to_num(flat_values, nan=no_data_value)
    else:
        raise ValueError(f"Property '{property_name}' has unsupported dtype {flat_values.dtype}")

    return {
        "values": np.asarray(flat_values, dtype=export_dtype),
        "storage_type": storage_type,
        "element_size": element_size,
        "no_data_value": no_data_value,
    }


def _write_structured_grid_gocad(grid, file_name):
    """Write a StructuredGrid to GOCAD VOXET format."""
    vo_path = Path(file_name).with_suffix(".vo")
    axis_n = np.asarray(grid.nsteps, dtype=int)
    property_source = grid.properties
    axis_min = np.min(grid.nodes, axis=0)
    axis_max = np.max(grid.nodes, axis=0)

    if property_source:
        if grid.cell_properties:
            logger.warning(
                "StructuredGrid GOCAD export uses point properties; cell_properties were not exported"
            )
    elif grid.cell_properties:
        axis_n = np.asarray(grid.nsteps, dtype=int) - 1
        if np.any(axis_n <= 0):
            raise ValueError("StructuredGrid cell_properties require at least two grid nodes per axis")
        property_source = grid.cell_properties
        axis_min = np.min(grid.cell_centres, axis=0)
        axis_max = np.max(grid.cell_centres, axis=0)
    else:
        raise ValueError("StructuredGrid has no properties to export to GOCAD")

    export_properties = []
    for index, (property_name, values) in enumerate(property_source.items(), start=1):
        export_info = _normalise_voxet_property(values, property_name, axis_n)
        safe_name = re.sub(r"[^0-9A-Za-z_-]+", "_", property_name).strip("_") or f"property_{index}"
        data_path = vo_path.with_name(f"{vo_path.stem}_{safe_name}@@")
        export_properties.append(
            {
                "index": index,
                "name": property_name,
                "data_path": data_path,
                **export_info,
            }
        )

    with open(vo_path, "w") as fp:
        fp.write(
            f"""GOCAD Voxet 1
HEADER {{
name: {grid.name}
}}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME Default
AXIS_NAME \"X\" \"Y\" \"Z\"
AXIS_UNIT \"m\" \"m\" \"m\"
ZPOSITIVE Elevation
END_ORIGINAL_COORDINATE_SYSTEM
AXIS_O 0.000000 0.000000 0.000000
AXIS_U 1.000000 0.000000 0.000000
AXIS_V 0.000000 1.000000 0.000000
AXIS_W 0.000000 0.000000 1.000000
AXIS_MIN {axis_min[0]} {axis_min[1]} {axis_min[2]}
AXIS_MAX {axis_max[0]} {axis_max[1]} {axis_max[2]}
AXIS_N {axis_n[0]} {axis_n[1]} {axis_n[2]}
AXIS_NAME \"X\" \"Y\" \"Z\"
AXIS_UNIT \"m\" \"m\" \"m\"
AXIS_TYPE even even even
"""
        )
        for export_property in export_properties:
            fp.write(
                f"""PROPERTY {export_property['index']} {export_property['name']}
PROPERTY_CLASS {export_property['index']} {export_property['name']}
PROP_UNIT {export_property['index']} {export_property['name']}
PROPERTY_CLASS_HEADER {export_property['index']} {export_property['name']} {{
}}
PROPERTY_SUBCLASS {export_property['index']} QUANTITY {export_property['storage_type']}
"""
            )
            if export_property["no_data_value"] is not None:
                fp.write(
                    f"PROP_NO_DATA_VALUE {export_property['index']} {export_property['no_data_value']}\n"
                )
            fp.write(
                f"""PROP_ETYPE {export_property['index']} IEEE
PROP_FORMAT {export_property['index']} RAW
PROP_ESIZE {export_property['index']} {export_property['element_size']}
PROP_OFFSET {export_property['index']} 0
PROP_FILE {export_property['index']} {export_property['data_path'].name}
"""
            )
        fp.write("END\n")

    for export_property in export_properties:
        with open(export_property["data_path"], "wb") as fp:
            export_property["values"].tofile(fp)

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
    properties_header = None
    if surf.properties:

        properties_header = f"""PROPERTIES {' '.join(list(surf.properties.keys()))}
NO_DATA_VALUES -99999
PROPERTY_CLASSES {' '.join(list(surf.properties.keys()))}
        """

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
{properties_header if properties_header else ""}
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
                fd.write(f"VRTX {v_idx:} {vert[0]} {vert[1]} {vert[2]}")
                if surf.properties:
                    for value in surf.properties.values():
                        fd.write(f" {value[idx]}")
                fd.write("\n")
                v_map[idx] = v_idx
                v_idx += 1
        for face in surf.triangles:
            if face[0] in v_map and face[1] in v_map and face[2] in v_map:
                fd.write(f"TRGL {v_map[face[0]]} {v_map[face[1]]} {v_map[face[2]]} \n")
        fd.write("END\n")
    return True


# def _write_pointset(points, file_name):
#     """
#     Write out a GOCAD VS file for a pointset

#     Parameters
#     ----------
#     points: SimpleNamespace()
#         Details of the points, as a SimpleNamespace() object. Fields are:
#             locations: locations, numpy ndarray with dtype = float64 & shape = (N,3)
#             vectors: vectors, numpy ndarray with dtype = float64 & shape = (N,3)
#             name: name of feature e.g. fault or supergroup, string

#     file_name: string
#         Desired filename

#     Returns
#     -------
#     True if successful

#     """
#     file_name = Path(file_name).with_suffix(".vs")
#     with open(f"{file_name}", "w") as fd:
#         fd.write(
#             f"""GOCAD VSet 1
# HEADER {{
# name: {points.name}
# }}
# GOCAD_ORIGINAL_COORDINATE_SYSTEM
# NAME Default
# PROJECTION Unknown