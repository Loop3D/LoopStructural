"""
2b. Refolded folds
===================


"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer, RotationAnglePlotter
from LoopStructural.datasets import load_laurent2016
import pandas as pd

# logging.getLogger().setLevel(logging.INFO)

# load in the data from the provided examples
data, bb = load_laurent2016()
# bb[1,2] = 10000

data.head()

newdata = pd.DataFrame(
    [[5923.504395, 4748.135254, 3588.621094, "s2", 1.0]],
    columns=["X", "Y", "Z", "feature_name", "val"],
)
data = pd.concat([data, newdata], sort=False)

rotation = [-69.11979675292969, 15.704944610595703, 6.00014591217041]


######################################################################
# Modelling S2
# ~~~~~~~~~~~~
#

model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)
s2 = model.create_and_add_fold_frame("s2", nelements=10000, buffer=0.5, solver="lu", damp=True)
viewer = LavaVuModelViewer(model)
viewer.add_scalar_field(s2[0], cmap="prism")
viewer.add_isosurface(s2[0], slices=[0, 1])
viewer.add_data(s2[0])
viewer.rotate(rotation)
viewer.display()


######################################################################
# Modelling S1
# ~~~~~~~~~~~~
#

s1 = model.create_and_add_folded_fold_frame(
    "s1", av_fold_axis=True, nelements=50000, buffer=0.3, solver="lu"
)


viewer = LavaVuModelViewer(model)
viewer.add_scalar_field(s1[0], cmap="prism")
viewer.rotate([-69.11979675292969, 15.704944610595703, 6.00014591217041])
viewer.display()

######################################################################
# S2/S1 S-Plots
# ~~~~~~~~~~~~~
#
s2_s1_splot = RotationAnglePlotter(s1)
s2_s1_splot.add_fold_limb_data()
s2_s1_splot.add_fold_limb_curve()
# fig, ax = plt.subplots(1,2,figsize=(10,5))
# x = np.linspace(s2[0].min(),s2[0].max(),1000)
# ax[0].plot(x,s1['fold'].fold_limb_rotation(x))
# ax[0].plot(s1['fold'].fold_limb_rotation.fold_frame_coordinate,s1['fold'].fold_limb_rotation.rotation_angle,'bo')
# ax[1].plot(s1['limb_svariogram'].lags,s1['limb_svariogram'].variogram,'bo')


######################################################################
# Modelling S0
# ~~~~~~~~~~~~
#

s0 = model.create_and_add_folded_foliation(
    "s0",
    av_fold_axis=True,
    nelements=50000,
    buffer=0.2,
    damp=True,
    solver="lu",
)

viewer = LavaVuModelViewer(model)
viewer.add_scalar_field(s0, cmap="tab20")
viewer.rotate([-69.11979675292969, 15.704944610595703, 6.00014591217041])
viewer.display()

######################################################################
# S1/S0 S-Plots
# ~~~~~~~~~~~~~
#
s1_s0_splot = RotationAnglePlotter(s0)
s1_s0_splot.add_fold_limb_data()
s1_s0_splot.add_fold_limb_curve()

# fig, ax = plt.subplots(1,2,figsize=(10,5))
# x = np.linspace(s1[0].min(),s1[0].max(),1000)
# ax[0].plot(x,s0['fold'].fold_limb_rotation(x))
# ax[0].plot(s0['fold'].fold_limb_rotation.fold_frame_coordinate,s0['fold'].fold_limb_rotation.rotation_angle,'bo')
# ax[1].plot(s0['limb_svariogram'].lags,s1['limb_svariogram'].variogram,'bo')

viewer = LavaVuModelViewer(model)
viewer.add_isosurface(s0, nslices=10, paint_with=s0, cmap="tab20")
# viewer.add_data(s0)
# viewer.add_fold(s0['fold'],locations=s0['support'].barycentre[::80])
viewer.rotate([-69.11979675292969, 15.704944610595703, 6.00014591217041])
viewer.display()
