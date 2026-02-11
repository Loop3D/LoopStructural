import numpy as np
import pandas as pd
from LoopStructural import InterpolatorBuilder, BoundingBox
from LoopStructural.visualisation import Loop3DView

x = np.linspace(0, 1, 10)
y = np.zeros(10)
z = np.ones(10)
v = np.zeros(10)

data = pd.concat(
    [
        pd.DataFrame({"X": x, "Y": y, "Z": z, "val": v}),
        pd.DataFrame({"X": x, "Y": y + 0.5, "Z": z, "val": v + 0.5}),
    ]
)

bounding_box = BoundingBox(np.array([-2, -2, -2]), np.array([2, 2, 2]))

interpolator = (
    InterpolatorBuilder(
        interpolatortype="FDI", nelements=1e4, bounding_box=bounding_box
    )
    .add_value_constraints(data.values)
    .setup_interpolator()
    .build()
)
mesh = bounding_box.structured_grid()
mesh.properties["val"] = interpolator.evaluate_value(mesh.nodes)
view = Loop3DView()
view.add_mesh(mesh.vtk().contour([0, 0.5]), opacity=0.4)
view.remove_scalar_bar()
view.add_points(data[["X", "Y", "Z"]].values, scalars=data["val"].values)
view.remove_scalar_bar()

view.view_yz()
view.show()