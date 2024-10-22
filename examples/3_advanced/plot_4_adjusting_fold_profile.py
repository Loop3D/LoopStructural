"""
Custom fold profile
------------------------
"""

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fold.fold_function import (
    FoldRotationType,
    get_fold_rotation_profile,
)
import LoopStructural

LoopStructural.setLogging("warning")


from LoopStructural.visualisation import Loop3DView, RotationAnglePlotter


from LoopStructural.datasets import load_laurent2016
import pandas as pd

data, bb = load_laurent2016()
data.head()
newdata = pd.DataFrame(
    [[5923.504395, 4748.135254, 3588.621094, "s2", 1.0]],
    columns=["X", "Y", "Z", "feature_name", "val"],
)
data = pd.concat([data, newdata], sort=False)
rotation = [-69.11979675292969, 15.704944610595703, 6.00014591217041]
model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)
s2 = model.create_and_add_fold_frame(
    "s2",
    nelements=10000,
    buffer=0.5,
)
s1 = model.create_and_add_folded_fold_frame(
    "s1",
    av_fold_axis=True,
    nelements=1e4,
    limb_profile_type=FoldRotationType.FOURIER_SERIES,
)

profile = s1[0].builder.fold_limb_rotation
profile.plot()
from ipywidgets import interact


@interact(c0=(-90, 90, 0.01), c1=(-90, 90, 0.01), c2=(-90, 90, 0.01), w=(0, 20, 0.01))
def fn(c0=profile.c0, c1=profile.c1, c2=profile.c2, w=profile.w):
    profile.c0 = c0
    profile.c1 = c1
    profile.c2 = c2
    profile.w = w
    profile.plot()
    # s1[0].scalar_field().vtk().contour(5).plot()
