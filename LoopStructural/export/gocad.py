import numpy as np


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
#     from pathlib import Path

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
