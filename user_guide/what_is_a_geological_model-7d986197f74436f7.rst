import numpy as np
import pandas as pd
from LoopStructural.visualisation import Loop3DView
x = np.linspace(0, 1, 10)
y = np.zeros(10)
z = np.ones(10)
v = np.zeros(10)

data = pd.DataFrame({"X": x, "Y": y, "Z": z, "val": v})
vector = pd.DataFrame(
    [[0.5, 0, 1.0, 0, 1, 0]], columns=["X", "Y", "Z", "nx", "ny", "nz"]
)

view = Loop3DView()
view.add_points(data[["X", "Y", "Z"]].values)
view.add_arrows(
    vector[["X", "Y", "Z"]].values,
    vector[["nx", "ny", "nz"]].values,
)
view.show()