"""
2b. Refolded folds
===================
The previous example modelled a single fold generation using a fold
frame. Multiply-deformed terranes often contain **refolded folds**, where
an earlier folded foliation is itself folded by a later deformation event.
LoopStructural handles this by nesting fold frames: a fold frame can
itself be folded by an older fold frame, and a folded foliation can then
be built within that nested coordinate system.

This example builds three progressively older/more-deformed features from
the Laurent et al. (2016) synthetic refolded-fold dataset:

* ``s2`` - the youngest fold frame, built directly from the data
* ``s1`` - an older fold frame, itself folded within ``s2``
* ``s0`` - the original bedding, folded within ``s1``
"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import Loop3DView, RotationAnglePlotter
from LoopStructural.datasets import load_laurent2016
import pandas as pd

data, bb = load_laurent2016()
data.head()

# add an extra value constraint for "s2" so that its scalar field has at
# least two distinct values to interpolate between
newdata = pd.DataFrame(
    [[5923.504395, 4748.135254, 3588.621094, "s2", 1.0]],
    columns=["X", "Y", "Z", "feature_name", "val"],
)
data = pd.concat([data, newdata], sort=False)

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data

######################################################################
# Modelling S2
# ~~~~~~~~~~~~
# ``s2`` is the youngest, least-deformed fold generation, so it can be
# built as a standard fold frame directly from the orientation and
# lineation observations.

s2 = model.create_and_add_fold_frame("s2", nelements=10000, buffer=0.5, solver="lu", damp=True)
viewer = Loop3DView(model)
viewer.plot_scalar_field(s2[0], cmap="prism")
viewer.plot_surface(s2[0], [0, 1])
viewer.plot_data(s2[0])
viewer.display()


######################################################################
# Modelling S1
# ~~~~~~~~~~~~
# ``s1`` is an older fold frame that has itself been refolded by the
# ``s2`` deformation event, so it is built with
# :code:`create_and_add_folded_fold_frame`, passing ``s2`` as the fold
# frame it is folded within, rather than the plain
# :code:`create_and_add_fold_frame` used for ``s2`` above.

s1 = model.create_and_add_folded_fold_frame(
    "s1", fold_frame=s2, av_fold_axis=True, nelements=50000, buffer=0.3, limb_wl=4
)

viewer = Loop3DView(model)
viewer.plot_scalar_field(s1[0], cmap="prism")
viewer.display()

######################################################################
# S2/S1 S-Plots
# ~~~~~~~~~~~~~
# The fold limb rotation angle of ``s1`` plotted against the ``s2`` fold
# frame coordinate - the same rotation-angle vs coordinate relationship
# used in the single-fold example, just calculated within the nested
# frame.
s2_s1_splot = RotationAnglePlotter(s1)
s2_s1_splot.add_fold_limb_data()
s2_s1_splot.add_fold_limb_curve()


######################################################################
# Modelling S0
# ~~~~~~~~~~~~
# ``s0`` is the original bedding, folded within the (already refolded)
# ``s1`` fold frame using :code:`create_and_add_folded_foliation`, in the
# same way the single-fold example folded ``s0`` within ``s1`` directly.

s0 = model.create_and_add_folded_foliation(
    "s0",
    fold_frame=s1,
    av_fold_axis=True,
    nelements=50000,
    buffer=0.2,
)

viewer = Loop3DView(model)
viewer.plot_scalar_field(s0, cmap="tab20")
viewer.display()

######################################################################
# S1/S0 S-Plots
# ~~~~~~~~~~~~~
s1_s0_splot = RotationAnglePlotter(s0)
s1_s0_splot.add_fold_limb_data()
s1_s0_splot.add_fold_limb_curve()

viewer = Loop3DView(model)
viewer.plot_surface(s0, 10, paint_with=s0, cmap="tab20")
viewer.display()
