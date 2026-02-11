vector[["nx", "ny", "nz"]] = vector[["nx", "ny", "nz"]]*2
bounding_box = BoundingBox(np.array([-2, -2, -2]), np.array([2, 2, 2]))

interpolator = (
    InterpolatorBuilder(
        interpolatortype="FDI", nelements=1e4, bounding_box=bounding_box
    ).add_value_constraints(data.values)
    .add_normal_constraints(vector.values).setup_interpolator()
    .build()
)
mesh2 = bounding_box.structured_grid()
mesh2.properties['val'] = interpolator.evaluate_value(mesh2.nodes)
view = Loop3DView()
view.add_mesh(mesh2.vtk().contour([-1, 0, 1]),color='red')
view.add_mesh(mesh.vtk().contour([-1, 0, 1]),color='blue')
view.add_points(data[["X", "Y", "Z"]].values)
view.add_arrows(
    vector[["X", "Y", "Z"]].values,
    vector[["nx", "ny", "nz"]].values,
)
view.remove_scalar_bar()
view.view_yz()
view.show()