"""
Building a model using the ProcessInputData
===========================================
There is a disconnect between the input data required by 3D modelling software and a geological map. In LoopStructural the geological model is a collection of implicit functions that can be mapped to the distribution of stratigraphic units and the location of fault surfaces. Each implicit function is approximated from the observations of the stratigraphy, this requires grouping conformable geological units together as a singla implicit function, mapping the different stratigraphic horizons to a value of the implicit function and determining the relationship with geological structures such as faults. In this tutorial the **ProcessInputData** class will be used to convert geologically meaningful datasets to input for LoopStructural. The **ProcessInputData** class uses:* stratigraphic contacts* stratigraphic orientations* stratigraphic thickness* stratigraphic order To build a model of stratigraphic horizons and:* fault locations* fault orientations * fault properties* fault edges To use incorporate faults into the geological model. """

##############################
# Imports
# ~~~~~~~


from LoopStructural.modelling import ProcessInputData, Map2LoopProcessor
from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##############################
# Read stratigraphy from csv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~


contacts = pd.read_csv("contacts.csv")
stratigraphic_orientations = pd.read_csv("stratigraphic_orientations.csv")
thicknesses = dict(
    zip(
        list(
            pd.read_csv(
                "stratigraphic_thickness.csv", skiprows=1, names=["name", "thickness"]
            )["name"]
        ),
        list(
            pd.read_csv(
                "stratigraphic_thickness.csv", skiprows=1, names=["name", "thickness"]
            )["thickness"]
        ),
    )
)
stratigraphic_order = pd.read_csv("stratigraphic_order.csv")

##############################
#  Stratigraphic Contacts
# ***********************


contacts

fig, ax = plt.subplots(1)
ax.scatter(
    contacts["X"], contacts["Y"], c=contacts["name"].astype("category").cat.codes
)
ax.set_title("Contact data")

##############################
# Stratigraphic orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stratigraphic orientations needs to have X, Y, Z and either azimuth and dip, dipdirection and dip, strike and dip (RH thumb rule) or the vector components of the normal vector (nx, ny, nz)

stratigraphic_orientations

##############################
# Stratigraphic thickness
# ~~~~~~~~~~~~~~~~~~~~~~~
# Stratigraphic thickness should be a dictionary containing the unit name (which should be in the contacts table) and the corresponding thickness of this unit.

thicknesses

##############################
# Bounding box
# ~~~~~~~~~~~~
# * Origin - bottom left corner of the model # * Maximum - top right hand corner of the model

bbox = pd.read_csv(
    "bbox.csv", index_col=0, header=None, names=["X", "Y", "Z"]
)  # ,'r').read().split('\n')
origin = bbox.loc["origin"].to_numpy()  # np.array(bbox[0].split(',')[1:],dtype=float)
maximum = bbox.loc["maximum"].to_numpy()  # np.array(bbox[1].split(',')[1:],dtype=float)

bbox

##############################
# Stratigraphic column
# ~~~~~~~~~~~~~~~~~~~~
# The order of stratrigraphic units is defined a list of tuples containing the name of the group and the order of units within the group. For example there are 7 units in the following example that form two groups.

# example nested list
[
    ("youngest_group", ["unit1", "unit2", "unit3", "unit4"]),
    ("older_group", ["unit5", "unit6", "unit7"]),
]

stratigraphic_order

order = [("supergroup_0", list(stratigraphic_order["unit name"]))]

##############################
# Building a stratigraphic model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A ProcessInputData onject can be built from these datasets using the argument names. A full list of possible arguments can be found in the documentation.


processor = ProcessInputData(
    contacts=contacts,
    contact_orientations=stratigraphic_orientations.rename(
        {"formation": "name"}, axis=1
    ),
    thicknesses=thicknesses,
    stratigraphic_order=order,
    origin=origin,
    maximum=maximum,
)

##############################
# The process input data can be used to directly build a geological model

model = GeologicalModel.from_processor(processor)
model.update()

##############################
# Or build directly from the dataframe and processor attributes.

model2 = GeologicalModel(processor.origin, processor.maximum)
model2.data = processor.data
model2.create_and_add_foliation("supergroup_0")
model2.update()

##############################
# Visualising model
# ~~~~~~~~~~~~~~~~~


view = LavaVuModelViewer(model)
view.add_model_surfaces()
view.rotation = [-37.965614318847656, 13.706363677978516, 3.110347032546997]
view.display()

##############################
# Adding faults
# ~~~~~~~~~~~~~


fault_locations = pd.read_csv("fault_locations.csv")
fault_orientations = pd.read_csv("fault_orientations.csv")

fault_orientations

fault_edges = []
with open("fault_edges.txt", "r") as f:
    for l in f.read().split("\n"):
        faults = l.split(",")
        if len(faults) == 2:
            fault_edges.append((faults[0], faults[1]))

fault_edges

fault_properties = pd.read_csv("fault_displacement.csv", index_col=0)

fault_properties

processor = ProcessInputData(
    contacts=contacts,
    contact_orientations=stratigraphic_orientations.rename(
        {"formation": "name"}, axis=1
    ),
    thicknesses=thicknesses,
    stratigraphic_order=order,
    origin=origin,
    maximum=maximum,
    fault_edges=fault_edges,
    fault_orientations=fault_orientations,
    fault_locations=fault_locations,
    fault_properties=fault_properties,
)

model = GeologicalModel.from_processor(processor)
model.update()

view = LavaVuModelViewer(model)
view.add_model_surfaces()
view.rotation = [-37.965614318847656, 13.706363677978516, 3.110347032546997]
view.display()
