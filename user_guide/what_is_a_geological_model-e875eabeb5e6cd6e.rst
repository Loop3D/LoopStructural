view = Loop3DView()
view.add_mesh(mesh.vtk().contour([0]))
view.add_points(data[["X", "Y", "Z"]].values)
view.add_arrows(
    vector[["X", "Y", "Z"]].values,
    vector[["nx", "ny", "nz"]].values,
)
view.show()