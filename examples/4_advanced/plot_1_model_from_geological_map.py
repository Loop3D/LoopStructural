"""
4a. Building a model using ProcessInputData
=============================================
There is a disconnect between the input data required by 3D modelling
software and a geological map. In LoopStructural the geological model is
a collection of implicit functions that can be mapped to the
distribution of stratigraphic units and the location of fault surfaces.
Building each implicit function from raw map observations requires
grouping conformable geological units together as a single implicit
function, mapping the different stratigraphic horizons to a value of
that implicit function, and determining the relationship with geological
structures such as faults.

The **ProcessInputData** class automates this conversion from
geologically meaningful datasets to LoopStructural input. It uses:

* stratigraphic contacts, orientations, thickness and order - to build a
  model of the stratigraphic horizons, and
* fault locations, orientations, properties and edges - to incorporate
  faults into the geological model.
"""

##############################
# Imports
# ~~~~~~~


import matplotlib.pyplot as plt

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_geological_map_data
from LoopStructural.modelling import ProcessInputData
from LoopStructural.visualisation import Loop3DView

##############################
# Read stratigraphy from csv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

(
    contacts,
    stratigraphic_orientations,
    stratigraphic_thickness,
    stratigraphic_order,
    bbox,
    fault_locations,
    fault_orientations,
    fault_properties,
    fault_edges,
) = load_geological_map_data()

thicknesses = dict(
    zip(
        list(stratigraphic_thickness["name"]),
        list(stratigraphic_thickness["thickness"]),
    )
)

##############################
#  Stratigraphic Contacts
# ***********************


contacts

fig, ax = plt.subplots(1)
ax.scatter(contacts["X"], contacts["Y"], c=contacts["name"].astype("category").cat.codes)
ax.set_title("Contact data")
plt.show()

##############################
# Stratigraphic orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stratigraphic orientations needs to have X, Y, Z and either azimuth and dip, dipdirection and dip, strike
# and dip (RH thumb rule) or the vector components of the normal vector (nx, ny, nz)

stratigraphic_orientations

##############################
# Stratigraphic thickness
# ~~~~~~~~~~~~~~~~~~~~~~~
# Stratigraphic thickness should be a dictionary containing the unit name (which should be in the contacts table)
# and the corresponding thickness of this unit.

thicknesses

##############################
# Bounding box
# ~~~~~~~~~~~~
# * Origin - bottom left corner of the model
# * Maximum - top right hand corner of the model

origin = bbox.loc["origin"].to_numpy()
maximum = bbox.loc["maximum"].to_numpy()

bbox

##############################
# Stratigraphic column
# ~~~~~~~~~~~~~~~~~~~~
# The order of stratigraphic units is defined as a list of tuples
# containing the name of the group and the order of units within the
# group, oldest last. For example, the following would describe 7 units
# forming two groups::
#
#    [
#        ("youngest_group", ["unit1", "unit2", "unit3", "unit4"]),
#        ("older_group", ["unit5", "unit6", "unit7"]),
#    ]
#
# Here all the units belong to a single group, "supergroup_0", since the
# dataset only contains one conformable sequence.

stratigraphic_order

order = [("supergroup_0", list(stratigraphic_order["unit name"]))]

##############################
# Building a stratigraphic model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A ProcessInputData object can be built from these datasets using the
# argument names. A full list of possible arguments can be found in the
# documentation.


processor = ProcessInputData(
    contacts=contacts,
    contact_orientations=stratigraphic_orientations.rename({"formation": "name"}, axis=1),
    thicknesses=thicknesses,
    stratigraphic_order=order,
    origin=origin,
    maximum=maximum,
)
processor.foliation_properties["supergroup_0"] = {"regularisation": 1.0}
##############################
# ``GeologicalModel.from_processor`` builds a geological model directly
# from the processor - grouping the units, mapping them to scalar field
# values and adding the foliation for you.

model = GeologicalModel.from_processor(processor)
model.update()

##############################
# The same result can also be reached by hand, using the processor's
# ``data`` dataframe (already in LoopStructural's X/Y/Z/feature_name/val
# form) directly with the usual ``GeologicalModel`` API.

model2 = GeologicalModel(processor.origin, processor.maximum)
model2.data = processor.data
model2.create_and_add_foliation("supergroup_0")
model2.update()

##############################
# Visualising model
# ~~~~~~~~~~~~~~~~~

view = Loop3DView(model)
view.plot_model_surfaces()
view.display()

##############################
# Adding faults
# ~~~~~~~~~~~~~
# Faults are added to ``ProcessInputData`` the same way as the
# stratigraphy: ``fault_locations``/``fault_orientations`` give the
# geometry (analogous to ``contacts``/``contact_orientations``),
# ``fault_properties`` gives per-fault parameters like displacement, and
# ``fault_edges`` declares which faults interact with each other (see the
# fault network example in :code:`3_fault` for how the interaction angle
# is used).

fault_orientations

fault_edges

fault_properties

processor = ProcessInputData(
    contacts=contacts,
    contact_orientations=stratigraphic_orientations.rename({"formation": "name"}, axis=1),
    thicknesses=thicknesses,
    stratigraphic_order=order,
    origin=origin,
    maximum=maximum,
    fault_edges=fault_edges,
    fault_orientations=fault_orientations,
    fault_locations=fault_locations,
    fault_properties=fault_properties,
)
processor.foliation_properties['supergroup_0']['regularisation'] = 1.0
model = GeologicalModel.from_processor(processor)
model.update()

view = Loop3DView(model)
view.plot_model_surfaces()
view.display()

##############################
# Visualise stratigraphic column
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from LoopStructural.visualisation import StratigraphicColumnView

scv = StratigraphicColumnView(model)
scv.plot()
