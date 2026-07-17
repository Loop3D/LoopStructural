"""
1d. Multiple groups
===================
The previous examples in this section built a model with a single
conformable series. Most geological models require more than one series -
for example where an unconformity separates two packages of rocks that were
deposited or intruded at different times and are not conformable with each
other. Each of these packages needs its own implicit function ("group"),
and the relationship between groups (unconformable, intrusive, etc.) needs
to be defined explicitly.

This example reuses the Claudius dataset and splits it into two groups
separated by an unconformity.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView


data, bb = load_claudius()
data = data.reset_index()

# Flip the sign of the scalar field/normals so the "strati2" group (added
# below) increases in the opposite direction to "strati".
data.loc[:, "val"] *= -1
data.loc[:, ["nx", "ny", "nz"]] *= -1

# Manually reassign a single data point to a second feature, "strati2", so
# that there are observations available to constrain it independently of
# "strati".
data.loc[792, "feature_name"] = "strati2"
data.loc[792, ["nx", "ny", "nz"]] = [0, 0, 1]
data.loc[792, "val"] = 0

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data

######################################################################
# Adding an unconformity
# ~~~~~~~~~~~~~~~~~~~~~~
# ``strati2`` is added first and marked as unconformable using
# :code:`model.add_unconformity(feature, value)`. This tells the model
# that everything below the given isovalue of ``strati2`` belongs to an
# older, separately-interpolated package - which is what allows ``strati``
# to be built afterwards without being affected by the ``strati2``
# observations.

strati2 = model.create_and_add_foliation(
    "strati2",
    interpolatortype="FDI",
    nelements=1e4,
)
uc = model.add_unconformity(strati2, 1)

strati = model.create_and_add_foliation(
    "strati",
    interpolatortype="FDI",
    nelements=1e4,
)

######################################################################
# Visualising both groups
# ~~~~~~~~~~~~~~~~~~~~~~~
# Each group is a separate scalar field, so isosurfaces for each are added
# to the viewer independently.

viewer = Loop3DView(model)
viewer.plot_surface(strati2, value=[2, 1.5, 1])
viewer.plot_surface(strati, value=[0, -60, -250, -330], paint_with=strati)
viewer.display()
