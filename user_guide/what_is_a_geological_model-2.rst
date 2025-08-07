import numpy as np
import pandas as pd
from LoopStructural import InterpolatorBuilder, BoundingBox
from LoopStructural.visualisation import Loop3DView

bounding_box = BoundingBox(np.array([-2, -2, -2]), np.array([2, 2, 2]))

interpolator = (
    InterpolatorBuilder(
        interpolatortype="FDI", nelements=1e4, bounding_box=bounding_box
    ).add_value_constraints(data.values)
    .add_normal_constraints(vector.values).setup_interpolator()
    .build()
)
mesh = bounding_box.structured_grid()
mesh.properties['val'] = interpolator.evaluate_value(mesh.nodes)
view = Loop3DView()
view.add_mesh(mesh.vtk())
view.show()