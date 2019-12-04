from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer

import pandas as pd

data, boundary_points = load_noddy_single_fold()
data.head()

mdata = pd.concat([data[:30],data[data['type']=='s1']])
model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
model.set_model_data(mdata)
fold_frame = model.create_and_add_fold_frame('s1',nelements=10000)
stratigraphy = model.create_and_add_folded_foliation('s0',
                                               fold_frame,
                                                nelements=10000,
                                               fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                               limb_wl=1)
viewer = LavaVuModelViewer(model,background="white")
viewer.add_isosurface(fold_frame[0],colour='blue',isovalue=0.4,alpha=0.5)
viewer.add_isosurface(fold_frame[1],colour='green',alpha=0.5)


viewer.add_data(stratigraphy['feature'])
viewer.add_isosurface(stratigraphy['feature'],
                     voxet=model.voxet())
viewer.interactive()