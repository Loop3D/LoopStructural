15
#import the Forward Modelling Engine modules - LoopStructural
from LoopStructural import GeologicalModel
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
# from LoopStructural.visualisation.rotation_angle_plotter import RotationAnglePlotter
# import other libraries
import geopandas
import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt


import pandas as pd
boundary_points = np.zeros((2,3))
boundary_points[0,0] = 0
boundary_points[0,1] = 0
boundary_points[0,2] = 5000
boundary_points[1,0] = 10000
boundary_points[1,1] = 7000
boundary_points[1,2] = 10000

data = pd.read_pickle('../notebooks/onefolddata.pkl')

d = data[np.logical_or(data['random']<0.5,data['type']=='s1')]

folded_model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
folded_model.set_model_data(d)
fold_frame = folded_model.create_and_add_fold_frame('s1')
folded_foliation = folded_model.create_and_add_folded_foliation('s0',
                                                   fold_frame,
                                                   fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                                   limb_wl=8000,
                                                                nelements=100000)
viewer = LavaVuModelViewer(background="white")

# determine the number of unique surfaces in the model from the input data and then calculate isosurfaces for this
# unique = np.unique(stratigraphy_builder.interpolator.get_control_points()[:,3])
# viewer.add_scalar_field(folded_model.bounding_box,(38,55,30),
#                       'box',
#                      paint_with=folded_foliation,
#                      cmap='prism')
#     viewer.add_vector_data(stratigraphy.get_interpolator().get_norm_constraints()[:,:3],stratigraphy.interpolator.get_norm_constraints()[:,:3],'data')
viewer.add_data(folded_foliation)
viewer.add_isosurface(folded_foliation,
                       colour='purple',
#                        nslices=10,
                       # paint_with=f1_frame.features[0]
                       )
viewer.lv.rotate([-57.657936096191406, -13.939384460449219, -6.758780479431152])
viewer.interactive()
