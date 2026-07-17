"""
3c. Defining the fault displacement function
============================================
By default LoopStructural displaces a faulted surface following a smooth,
symmetric profile: displacement is greatest at the fault surface/centre
and decays to zero at the edges of the fault's ellipsoidal region of
influence, the same on both the hanging wall and footwall. Real faults
are not always this symmetric - for example, a **drag fault** shows extra
deformation of the faulted surface close to the fault on one side only.

This example shows how the default displacement profile looks (as three
1D functions of the three fault frame coordinates), then defines a custom
profile and uses it to build a drag fault.
"""

import numpy as np
import pandas as pd
import LoopStructural as LS

# A minimal dataset for a single vertical fault (two points defining its
# plane, "coord" 0 and 1) offsetting a single stratigraphic contact.

origin = [0, 0, 0]
extent = [10, 10, 10]

data = pd.DataFrame(
    [
        [5, 5, 5, 0, 0.70710678, 0.0, 0.70710678, 0, "fault"],
        [5, 5, 5, 0, -0.70710678, 0.0, 0.70710678, 1, "fault"],
        [8, 5, 5, 0, 0, 0, 1, np.nan, "strati"],
    ],
    columns=["X", "Y", "Z", "val", "nx", "ny", "nz", "coord", "feature_name"],
)

data

######################################################################
# Create model using the standard fault displacement model

model = LS.GeologicalModel(origin, extent)
model.data = data
model.create_and_add_fault(
    "fault",
    1,
    nelements=1000,
    interpolator_type="PLI",
    buffer=0.5,
    major_axis=10,
    minor_axis=3,
    intermediate_axis=10,
)
model.create_and_add_foliation(
    "strati", nelements=1000, interpolator_type="PLI", faults=[model["fault"]]
)


import LoopStructural.visualisation as vis

view = vis.Loop3DView(model)
view.plot_surface(model.features[0], value=[0])
view.plot_surface(model.features[1], value=5, paint_with=model.features[1])

view.display()

######################################################################
# The default displacement profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ``model['fault'].faultfunction`` holds three 1D profile functions, one
# per fault frame coordinate:
#
# * ``gx`` - displacement as a function of distance from the fault
#   surface (0 at the fault, decaying to 0 again at the edge of its
#   ellipsoidal region of influence). This is where the hanging
#   wall/footwall asymmetry lives - by default this profile is
#   antisymmetric, giving equal and opposite displacement on each side.
# * ``gy`` - displacement as a function of position along the fault slip
#   direction
# * ``gz`` - displacement as a function of position along the fault
#   extent (strike) direction
#
# The final displacement at a point is the product of all three profiles,
# scaled by the requested displacement magnitude. ``FaultDisplacement``
# has a convenience ``.plot()`` that shows all three together.

model['fault'].faultfunction.plot()

######################################################################
# A custom drag-fault profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A drag fault shows extra bending of the faulted surface adjacent to the
# fault on one side only. To reproduce this we replace the default,
# symmetric ``gx`` profile with a :code:`Composite` of a footwall function
# that decays away from the fault (as before) and a hanging wall function
# that is constant (:code:`Ones`) - i.e. no drag on the hanging wall side,
# full drag on the footwall side.

from LoopStructural.modelling.features.fault._fault_function import (
    FaultDisplacement,
    CubicFunction,
    Ones,
)

fw = CubicFunction()
fw.add_cstr(0, -1)
fw.add_grad(0, 0)
fw.add_cstr(-1, 0)
fw.add_grad(-1, 0)
fw.add_min(-1)
hw = Ones()
drag_fault = FaultDisplacement(hw=hw, fw=fw)

drag_fault.plot()

######################################################################
# Rebuilding the model with the custom profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The custom profile is passed in as ``faultfunction`` when the fault is
# created. The rest of the workflow is unchanged - "strati" still just
# needs to know which faults affect it.

model = LS.GeologicalModel(origin, extent)
model.data = data
model.create_and_add_fault(
    "fault",
    -1,
    nelements=1000,
    interpolator_type="PLI",
    buffer=0.5,
    major_axis=10,
    minor_axis=6,
    intermediate_axis=10,
    faultfunction=drag_fault,
)
model.create_and_add_foliation(
    "strati", nelements=1000, interpolator_type="PLI", faults=[model["fault"]]
)

view = vis.Loop3DView(model)
model.bounding_box.nelements = 1e5
view.plot_surface(model.features[0], value=[0])
view.plot_surface(model['strati'], value=5)

view.display()
