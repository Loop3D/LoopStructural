"""
1f. Radial Basis Interpolation Using Surfe
==========================================
Surfe is a C++ library written by Michael Hillier from the Geological Survey of Canada.
The library can be found on `<https://github.com/MichaelHillier/surfe>`_

Python wrappers allow for the library to be used from within a python interface. 
There are two ways of installing surfe:
1. Follow the CMake instructions on the surfe repository to build the code from source. 
2. Download the precompiled library from the github releases `<https://github.com/MichaelHillier/surfe/releases>`_

Once you have the correct file for your system you will need to add an environment variable
to windows :code:`SURFE`

"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer

from LoopStructural.datasets import load_claudius  # demo data

import numpy as np
import logging

logging.getLogger().setLevel(logging.INFO)

data, bb = load_claudius()  # claudius.get_data()
bb[1, 0] += 200
bb[0, 0] -= 200
bb[1, 1] += 200
bb[0, 1] -= 200
bb[1, 2] += 200
bb[0, 2] -= 200

model = GeologicalModel(bb[0, :], bb[1, :])
data["random"] = np.random.random(data.shape[0])
model.set_model_data(data[data["random"] < 0.01])  # [np.isnan(data['val'])])
strati = model.create_and_add_foliation("strati", interpolatortype="surfe", method="single_surface")
print(strati.evaluate_value(model.regular_grid((10, 10, 10))))
viewer = LavaVuModelViewer(model, background="white")

# determine the number of unique surfaces in the model from
# the input data and then calculate isosurfaces for this
unique = np.unique(strati.interpolator.get_value_constraints()[:, 3])
viewer.add_isosurface(model.features[0], slices=unique, cmap="prism", paint_with=model.features[0])
#
# # viewer.add_section(model.features[0],
# #                    axis='x',
# #                    value=0,
# #                    boundary_points=model.bounding_box,
# #                    nsteps=np.array([30,30,30]),
# #                     ,
# #                   cmap='prism')
# viewer.add_scalar_field(model.features[0],
#                      cmap='prism')

# # Add the data addgrad/addvalue arguments are optional
# viewer.add_data(model.features[0],addgrad=True,addvalue=True, cmap='prism')
# viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.interactive()  # to add an interactive display
