"""
3a. Modelling faults using structural frames
========================================

"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import Loop3DView
from LoopStructural.datasets import load_intrusion
import numpy as np
import matplotlib.pyplot as plt

data, bb = load_intrusion()


######################################################################
# Modelling faults using structural frames
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Standard implicit modelling techniques either treat faults as domain
# boundaries or use a step function in the implicit function to capture
# the displacement in the faulted surface.
#
# Adding faults into the implicit function using step functions is limited
# because this does not capture the kinematics of the fault. It
# effectively defines the fault displacement by adding a value to the
# scalar field on the hanging wall of the fault. In the example below a
# 2-D ellipsoidal function is combined with a step function to show how
# the resulting geometry results in a shrinking shape. This would be
# representative of modelling an intrusion.
#

intrusion = lambda x, y: (x * 2) ** 2 + (y**2)
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
xx, yy = np.meshgrid(x, y)
fault = np.zeros(xx.shape)
fault[yy > 0] = 50
val = intrusion(xx, yy) + fault


plt.contourf(val)


######################################################################
# LoopStructural applies structural frames to the fault geometry to
# capture the geometry and kinematics of the fault. A fault frame
# consisting of the fault surface, fault slip direction and fault extent
# are built from observations. The geometry of the deformed surface is
# then interpolated by first restoring the observations by combining the
# fault frame and an expected displacement model.
#

model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)
fault = model.create_and_add_fault(
    "fault", 500
)

viewer = Loop3DView(model)
viewer.plot_surface(
    fault,
    value=0,
    #                       slices=[0,1]#nslices=10
)
xyz = model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].to_numpy()
xyz = xyz[fault.evaluate(xyz).astype(bool), :]
viewer.plot_vector_field(fault)
viewer.add_points(
    model.rescale(
        model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values,
        inplace=False,
    ),
    name="prefault",
)
viewer.display()


displacement = 400  # INSERT YOUR DISPLACEMENT NUMBER HERE BEFORE #

model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)
fault = model.create_and_add_fault(
    "fault", displacement, nelements=2000, 
)
strati = model.create_and_add_foliation("strati")
model.update()
viewer = Loop3DView(model)
viewer.plot_surface(strati, value=0.)
# viewer.add_data(model.features[0][0])
viewer.plot_data(strati)
viewer.plot_surface(
    fault,
    value=0.,
    #                       slices=[0,1]#nslices=10
)
viewer.add_points(
    model.rescale(
        model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values,
        inplace=False,
    ),
    name="prefault",
)
viewer.display()
