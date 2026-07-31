"""
3a. Modelling faults using structural frames
=============================================
This tutorial introduces how LoopStructural represents faults, and
compares that to the simpler step-function approach used by many implicit
modelling tools.
"""

import matplotlib.pyplot as plt
import numpy as np

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_intrusion
from LoopStructural.visualisation import Loop3DView

data, bb = load_intrusion()


######################################################################
# Why not just use a step function?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Standard implicit modelling techniques either treat faults as domain
# boundaries or use a step function in the implicit function to capture
# the displacement of the faulted surface.
#
# Adding faults into the implicit function using step functions is limited
# because this does not capture the kinematics of the fault. It
# effectively defines the fault displacement by adding a value to the
# scalar field on the hanging wall of the fault. In the example below a
# 2-D ellipsoidal function is combined with a step function to show how
# the resulting geometry results in a shrinking shape rather than a
# displaced one - a step function on its own cannot reproduce a fault
# that both offsets *and* preserves the shape of a surface, which is what
# real faults do.

intrusion = lambda x, y: (x * 2) ** 2 + (y**2)
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
xx, yy = np.meshgrid(x, y)
fault = np.zeros(xx.shape)
fault[yy > 0] = 50
val = intrusion(xx, yy) + fault

plt.contourf(val)
plt.title("Step function added to an ellipsoidal field - shrinks, doesn't displace")
plt.show()


######################################################################
# Faults as structural frames
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LoopStructural instead applies structural frames to the fault geometry to
# capture the geometry and kinematics of the fault. A fault frame
# consisting of the fault surface, fault slip direction and fault extent
# are built from observations. The geometry of the deformed surface is
# then interpolated by first restoring the observations by combining the
# fault frame and an expected displacement model - i.e. undoing the fault
# to interpolate the surface, then reapplying the displacement.
#
# ``create_and_add_fault(name, displacement)`` is all that's needed to
# add a fault - ``displacement`` sets the maximum offset across the fault.

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
fault = model.create_and_add_fault("fault", 500)

viewer = Loop3DView(model)
viewer.plot_surface(fault, value=0)
viewer.plot_vector_field(fault)
viewer.add_points(
    model.rescale(
        model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values,
        inplace=False,
    ),
    name="prefault",
)
viewer.display()

######################################################################
# Faulting a stratigraphic surface
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding a foliation *after* the fault (as in the previous 1_basic
# examples) means it is automatically restored/displaced using the fault
# frame built above. Try changing ``displacement`` below and re-running to
# see how the offset of "strati" across the fault surface changes.

displacement = 400

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
fault = model.create_and_add_fault("fault", displacement, nelements=2000)
strati = model.create_and_add_foliation("strati")
model.update()

viewer = Loop3DView(model)
viewer.plot_surface(strati, value=0.0)
viewer.plot_data(strati)
viewer.plot_surface(fault, value=0.0)
viewer.add_points(
    model.rescale(
        model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values,
        inplace=False,
    ),
    name="prefault",
)
viewer.display()
