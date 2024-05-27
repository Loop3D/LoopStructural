"""
============================
1f. Unconformities
============================
This tutorial will demonstrate how to add unconformities to a mode using LoopStructural.

"""

from LoopStructural import GeologicalModel
import pandas as pd
import numpy as np
from LoopStructural.visualisation import Loop3DView

##################################################################################################
# Generate synthetic data
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Model 3 scalar fields where the top is horizontal, the middle is dipping and the bottom is horizontal.
data = pd.DataFrame(
    [
        [0, 0, 3, 0, 0, 0, 1, "unit_a"],
        [0, 0, 1, 1, 0, 0.7, 0.7, "unit_b"],
        [0, 0, 0, 0, 0, 0.7, 0.7, "unit_b"],
        [0, 0, -3, 0, 0, 0, -1, "unit_c"],
    ],
    columns=["X", "Y", "Z", "val", "nx", "ny", "nz", "feature_name"],
)

model = GeologicalModel(np.ones(3) * -5, np.ones(3) * 7)
model.data = data
model.create_and_add_foliation("unit_a")
model.create_and_add_foliation("unit_b")
model.create_and_add_foliation("unit_c")

model.update()
##################################################################################################
# Visualise the model without unconformities
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

view = Loop3DView(model)
view.plot_surface(model["unit_a"], value=5)
view.plot_surface(model["unit_b"], value=5)
view.plot_surface(model["unit_c"], value=5)
view.display()

##################################################################################################
# Add unconformities
# ~~~~~~~~~~~~~~~~~~
# We add two unconformities to the model
# 1. the isovalue of 0 of unit_a is an unconformity
# 2. the isovalue of 0 of unit_b is an unconformity
#
# This means unit_a should not occur below isovalue of 0,
# unit_b should truncate at unit_a isovalue 0 and
# unit_b should not occur below isovalue of 0
# and unit_c should not occur below unit_b isovalue of 0

model = GeologicalModel(np.ones(3) * -5, np.ones(3) * 7)
model.data = data
model.create_and_add_foliation("unit_a")
model.add_unconformity(model["unit_a"], 0)
model.create_and_add_foliation("unit_b")
model.add_unconformity(model["unit_b"], 0)
model.create_and_add_foliation("unit_c")

##################################################################################################
# We can examine the model by printing the object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(model)

model.update()

##################################################################################################
# Visualise the model without unconformities
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

view = Loop3DView(model)
view.plot_surface(model["unit_a"], value=5)
view.plot_surface(model["unit_b"], value=5)
view.plot_surface(model["unit_c"], value=5)
view.display()


##################################################################################################
# Adding onlap unconformity
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# We can also add onlap unconformities to the model, using the previous example lets change the unconformity
# between b and c to be an onlap. This means the geometry of c truncates b


model = GeologicalModel(np.ones(3) * -5, np.ones(3) * 7)
model.data = data
model.create_and_add_foliation("unit_a")
model.add_unconformity(model["unit_a"], 0)
model.create_and_add_foliation("unit_b")
model.create_and_add_foliation("unit_c")
model.add_onlap_unconformity(model["unit_c"], 0)

model.update()

##################################################################################################
# Visualise the model with onlap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

view = Loop3DView(model)
view.plot_surface(model["unit_a"], value=5)
view.plot_surface(model["unit_b"], value=5)
view.plot_surface(model["unit_c"], value=5)


view.display()
