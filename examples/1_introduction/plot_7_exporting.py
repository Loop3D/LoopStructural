"""

Exporting models 
===============================

Models can be exported to vtk, gocad and geoh5 formats.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

data, bb = load_claudius()

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation("strati")


######################################################################
# Export surfaces to vtk
# ~~~~~~~~~~~~~~~~~~~~~~
# Isosurfaces can be extracted from a geological feature by calling
# the `.surfaces` method on the feature. The argument for this method
# is the value, values or number of surfaces that are extracted.
# This returns a list of `LoopStructural.datatypes.Surface` objects
# These objects can be interrogated to return the triangles, vertices
# and normals. Or can be exported into another format using the `save`
# method. The supported file formats are `vtk`, `ts` and `geoh5`.
#

surfaces = model['strati'].surfaces(value=0.0)

print(surfaces)

print(surfaces[0].vtk)

# surfaces[0].save('text.geoh5')

######################################################################
# Export the model to geoh5
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# The entire model can be exported to a geoh5 file using the `save_model`
# method. This will save all the data, foliations, faults and other objects
# in the model to a geoh5 file. This file can be loaded into LoopStructural

# model.save('model.geoh5')
