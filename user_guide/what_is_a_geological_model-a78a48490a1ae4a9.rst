view = Loop3DView()
view.add_mesh(mesh.vtk().contour([-1, 0, 1]))
view.remove_scalar_bar()
view.add_points(data[["X", "Y", "Z"]].values)
view.add_arrows(
    vector[["X", "Y", "Z"]].values,
    vector[["nx", "ny", "nz"]].values,
)
view.remove_scalar_bar()
view.view_yz()
view.show()