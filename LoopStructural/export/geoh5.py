import geoh5py
import geoh5py.workspace
import numpy as np
from LoopStructural.datatypes import ValuePoints, VectorPoints, Surface, StructuredGrid
from typing import Union

def add_surface_to_geoh5(filename:str, surface:Surface, overwrite:bool=True, groupname:str="Loop"):
    """Add a surface to a geoh5 file.

    Parameters
    ----------
    filename : _type_
        path to the geoh5 file
    surface : Surface
        the surface to add to the geoh5 file
    overwrite : bool, optional
        whether to overwrite the object, by default True
    groupname : str, optional
        Geoscience analyst group name, by default "Loop"
    """
    with geoh5py.workspace.Workspace(filename) as workspace:
        group = workspace.get_entity(groupname)[0]
        if not group:
            group = geoh5py.groups.ContainerGroup.create(
                workspace, name=groupname, allow_delete=True
            )
        if surface.name in workspace.list_entities_name.values():
            existing_surf = workspace.get_entity(surface.name)
            existing_surf[0].allow_delete = True
            if overwrite:
                workspace.remove_entity(existing_surf[0])
        data = {}
        if surface.properties is not None:
            for k, v in surface.properties.items():
                data[k] = {'association': "VERTEX", "values": v}
        surface = geoh5py.objects.Surface.create(
            workspace,
            name=surface.name,
            vertices=surface.vertices,
            cells=surface.triangles,
            parent=group,
        )
        surface.add_data(data)


def add_points_to_geoh5(filename:str, point:Union[ValuePoints, VectorPoints], overwrite:bool=True, groupname:str="Loop") -> None:
    """Add points to a geoh5 file.
    Parameters
    ----------
    filename : str
        path to the geoh5 file
    point : ValuePoints or VectorPoints
        the points to add to the geoh5 file
    overwrite : bool, optional
        whether to overwrite the object, by default True
    groupname : str, optional
        Geoscience analyst group name, by default "Loop"
    """
    with geoh5py.workspace.Workspace(filename) as workspace:
        group = workspace.get_entity(groupname)[0]
        if not group:
            group = geoh5py.groups.ContainerGroup.create(
                workspace, name=groupname, allow_delete=True
            )
        if point.name in workspace.list_entities_name.values():
            existing_point = workspace.get_entity(point.name)
            existing_point[0].allow_delete = True
            if overwrite:
                workspace.remove_entity(existing_point[0])
        data = {}
        if point.properties is not None:
            for k, v in point.properties.items():
                try:
                    v = np.array(v)
                except ValueError:
                    raise ValueError(
                        f"Cannot convert {k} to numpy array. Please check the data type."
                    )
                if isinstance(v, np.ndarray):
                    v = v.flatten()
                
                data[k] = {'association': "VERTEX", "values": v}
        if isinstance(point, VectorPoints):
            if point.vectors is None:
                raise ValueError("Vectors cannot be None for VectorPoints.")
            point.vectors = np.array(point.vectors)
            if point.vectors.shape[1] != 3:
                raise ValueError(
                    f"Vectors should have shape (n_points, 3), but got {point.vectors.shape}."
                )
            if point.vectors.shape[0] != point.locations.shape[0]:
                raise ValueError(
                    f"Vectors and locations should have the same number of points, but got {point.vectors.shape[0]} and {point.locations.shape[0]}."
                )
            
            if isinstance(point.vectors, np.ndarray):
                point.vectors = point.vectors.astype(np.float64)

            if point.vectors.dtype != np.float64:
                raise ValueError(
                    f"Vectors should be of type float64, but got {point.vectors.dtype}."
                )
            
            data['vx'] = {'association': "VERTEX", "values": point.vectors[:, 0]}
            data['vy'] = {'association': "VERTEX", "values": point.vectors[:, 1]}
            data['vz'] = {'association': "VERTEX", "values": point.vectors[:, 2]}

        if isinstance(point, ValuePoints):
            data['val'] = {'association': "VERTEX", "values": point.values}
        point = geoh5py.objects.Points.create(
            workspace,
            name=point.name,
            vertices=point.locations,
            parent=group,
        )
        point.add_data(data)


def add_structured_grid_to_geoh5(filename:str, structured_grid:StructuredGrid, overwrite:bool=True, groupname:str="Loop"):
    """Add a structured grid to a geoh5 file.
    Parameters
    ----------
    filename : str
        path to the geoh5 file
    structured_grid : StructuredGrid
        the structured grid to add to the geoh5 file
    overwrite : bool, optional
        whether to overwrite the object, by default True
    groupname : str, optional
        Geoscience analyst group name, by default "Loop"
    """
    with geoh5py.workspace.Workspace(filename) as workspace:
        group = workspace.get_entity(groupname)[0]
        if not group:
            group = geoh5py.groups.ContainerGroup.create(
                workspace, name=groupname, allow_delete=True
            )
        if structured_grid.name in workspace.list_entities_name.values():
            existing_block = workspace.get_entity(structured_grid.name)
            existing_block[0].allow_delete = True
            if overwrite:
                workspace.remove_entity(existing_block[0])
        data = {}
        if structured_grid.cell_properties is not None:
            for k, v in structured_grid.cell_properties.items():
                data[k] = {
                    "association": "CELL",
                    "values": np.flipud(
                        np.rot90(v.reshape(structured_grid.nsteps - 1, order="F"), 1)
                    ).flatten(),
                }
        block = geoh5py.objects.BlockModel.create(
            workspace,
            name=structured_grid.name,
            origin=structured_grid.origin,
            u_cell_delimiters=np.cumsum(
                np.ones(structured_grid.nsteps[0]) * structured_grid.step_vector[0]
            ),  # Offsets along u
            v_cell_delimiters=np.cumsum(
                np.ones(structured_grid.nsteps[1]) * structured_grid.step_vector[1]
            ),  # Offsets along v
            z_cell_delimiters=np.cumsum(
                np.ones(structured_grid.nsteps[2]) * structured_grid.step_vector[2]
            ),  # Offsets along z (down)
            rotation=0.0,
            parent=group,
        )
        block.add_data(data)
