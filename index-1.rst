from LoopStructural import GeologicalModel
from LoopStructural.datatypes import BoundingBox
from LoopStructural.visualisation import Loop3DView
from LoopStructural.datasets import load_claudius

import numpy as np
data, bb = load_claudius()

#bb constaints origin and maximum of axis aligned bounding box
#data is a pandas dataframe with X,Y,Z,val,nx,ny,nz, feature_name

model = GeologicalModel(bb[0,:],bb[1,:])
model.data = data
# nelements specifies the number of discrete interpolation elements
# 'stratí' is the feature name in the data dataframe
model.create_and_add_foliation('strati',nelements=1e5)
model.update()
# get the value of the interpolator at some random locations
locations = np.array(
    [
        np.random.uniform(bb[0, 0], bb[1, 0],5),
        np.random.uniform(bb[0, 1], bb[1, 1],5),
        np.random.uniform(bb[0, 2], bb[1, 2],5),
    ]
).T
val = model.evaluate_feature_value('strati', locations)
# get the gradient of the interpolator
gradient = model.evaluate_feature_gradient('strati',locations)

#Plot the scalar field of the model
model['strati'].scalar_field().plot()