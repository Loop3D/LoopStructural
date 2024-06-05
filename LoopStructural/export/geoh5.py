import geoh5py
import numpy as np


def add_surface_to_geoh5(filename, surface, overwrite=True, groupname="Loop"):
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


def add_points_to_geoh5(filename, points, overwrite=True, groupname="Loop"):
    with geoh5py.workspace.Workspace(filename) as workspace:
        group = workspace.get_entity(groupname)[0]
        if not group:
            group = geoh5py.groups.ContainerGroup.create(
                workspace, name=groupname, allow_delete=True
            )
        for point in points:
            if point.name in workspace.list_entities_name.values():
                existing_point = workspace.get_entity(point.name)
                existing_point[0].allow_delete = True
                if overwrite:
                    workspace.remove_entity(existing_point[0])
            data = {}
            if point.properties is not None:
                for k, v in point.properties.items():
                    data[k] = {'association': "VERTEX", "values": v}
            point = geoh5py.objects.Points.create(
                workspace,
                name=point.name,
                vertices=point.vertices,
                parent=group,
            )
            point.add_data(data)


def add_block_model_to_geoh5py(filename, block_model, overwrite=True, groupname="Loop"):
    with geoh5py.workspace.Workspace(filename) as workspace:
        group = workspace.get_entity(groupname)[0]
        if not group:
            group = geoh5py.groups.ContainerGroup.create(
                workspace, name=groupname, allow_delete=True
            )
        if block_model.name in workspace.list_entities_name.values():
            existing_block = workspace.get_entity(block_model.name)
            existing_block[0].allow_delete = True
            if overwrite:
                workspace.remove_entity(existing_block[0])
        data = {}
        if block_model.properties is not None:
            for k, v in block_model.properties.items():
                data[k] = {'association': "CELL", "values": v}
        block = geoh5py.objects.BlockModel.create(
            workspace,
            name=block_model.name,
            origin=[25, -100, 50],
            u_cell_delimiters=np.cumsum(np.ones(16) * 5),  # Offsets along u
            v_cell_delimiters=np.cumsum(np.ones(32) * 5),  # Offsets along v
            z_cell_delimiters=np.cumsum(np.ones(16) * -2.5),  # Offsets along z (down)
            rotation=30.0,
            parent=group,
        )
        block.add_data(data)
