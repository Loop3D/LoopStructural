from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer

from LoopStructural.datasets import load_claudius #demo data

import pandas as pd
import glob
import numpy as np
import logging
logging.getLogger().setLevel(logging.INFO)

data, bb = load_claudius()#claudius.get_data()
bb[1,0]+=200
bb[0,0]-=200
bb[1,1]+=200
bb[0,1]-=200
bb[1,2]+=200
bb[0,2]-=200

model = GeologicalModel(bb[0,:],bb[1,:])
data['random'] = np.random.random(data.shape[0])
model.set_model_data(data[data['random'] < 0.01])#[np.isnan(data['val'])])
strati = model.create_and_add_foliation("strati",
                                           interpolatortype="surfe",
                                        )
viewer = LavaVuModelViewer(model,background="white")

# determine the number of unique surfaces in the model from
# the input data and then calculate isosurfaces for this
unique = np.unique(strati['feature'].interpolator.get_value_constraints()[:,3])
viewer.add_isosurface(model.features[0],
                       slices=unique,
                       cmap='prism',
                                         voxet=model.voxet(),
                      paint_with=model.features[0])

# viewer.add_section(model.features[0],
#                    axis='x',
#                    value=0,
#                    boundary_points=model.bounding_box,
#                    nsteps=np.array([30,30,30]),
#                    voxet=model.voxet(),
#                   cmap='prism')
viewer.add_scalar_field(model.features[0],
                     cmap='prism')

# Add the data addgrad/addvalue arguments are optional
viewer.add_data(model.features[0],addgrad=True,addvalue=True, cmap='prism')
viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.interactive()# to add an interactive display
