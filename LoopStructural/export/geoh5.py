import geoh5py


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


# def add_points_to_geoh5(filename, points, overwrite=True, groupname="Loop"):
