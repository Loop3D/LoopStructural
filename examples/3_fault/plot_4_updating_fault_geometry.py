"""
3d. Updating fault geometry
============================
Building a model can be expensive, so LoopStructural avoids
re-interpolating a feature until it's actually needed. Changing a
parameter on a feature's builder - for example a fault's ``minor_axis``,
the size of its ellipsoidal region of influence - just marks that feature
(and anything downstream of it, like a faulted foliation) as out of date;
the next call to :code:`model.update()`, or the next time the feature is
evaluated, transparently triggers a rebuild using the new parameter. You
don't need to rebuild the ``GeologicalModel`` from scratch to try out a
different parameter value.

This example builds a faulted model once, then changes the fault's
``minor_axis`` and rebuilds, comparing the two results.
"""

import numpy as np
import pandas as pd

import LoopStructural as LS
import LoopStructural.visualisation as vis

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

# The prepared example data is used below to build the faulted model.

######################################################################
# Build the model once

model = LS.GeologicalModel(origin, extent)
model.data = data
model.create_and_add_fault(
    "fault",
    10,
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
model.update()

point = np.array([[6.5, 5, 5]])
print(f"minor_axis={model['fault'].builder.fault_minor_axis}, "
      f"strati value at {point[0]}: {model['strati'].evaluate_value(point)[0]:.3f}")

view = vis.Loop3DView(model)
view.plot_surface(model['fault'], value=[0])
view.plot_surface(model['strati'], value=5, paint_with=model['strati'])
view.display()

######################################################################
# Change a fault parameter and update
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setting ``fault_minor_axis`` directly on the builder is all that's
# needed - :code:`model.update()` picks up the change and only
# re-interpolates the fault and the features that depend on it, rather
# than the whole model.

model['fault'].builder.fault_minor_axis = 6.0
model.update()

print(f"minor_axis={model['fault'].builder.fault_minor_axis}, "
      f"strati value at {point[0]}: {model['strati'].evaluate_value(point)[0]:.3f}")

view = vis.Loop3DView(model)
view.plot_surface(model['fault'], value=[0])
view.plot_surface(model['strati'], value=5, paint_with=model['strati'])
view.display()
