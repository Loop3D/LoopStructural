from skimage.measure import marching_cubes
import pandas as pd
import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

def calculate_fault_intersections(model):
    fault_names = []
    for f in model.features:
        if f.type == 'fault':
            fault_names.append(f.name)
    fault_matrix = pd.DataFrame(columns=fault_names)
    for f in fault_names:
        fault_matrix.loc[f,:] = 0
        
    for name in fault_names:
        xyz=model.regular_grid(shuffle=False)
        vals = model[name].evaluate_value(xyz)
        model.step_vector = (model.bounding_box[1,:]-model.bounding_box[0,:])/model.nsteps
        try:
                
            verts, faces, normals, values = marching_cubes(
                            vals.reshape(model.nsteps, order='C'),
                            0,
                            spacing=model.step_vector)
            verts = model.rescale(verts)
        except (ValueError, RuntimeError) as e:
                print(e)
                logger.warning("Cannot isosurface {} at {}, skipping".format(name,0))
                continue
        for name2 in fault_names:
            if name2  == name:
                continue
            val = model[name2].evaluate_value(model.scale(verts,inplace=False))

            if np.all(np.isnan(val)):
                continue
            if np.all(val[~np.isnan(val)] > 0):
                continue#print('all > 0')
            elif np.all(val[~np.isnan(val)] < 0):
                continue
            else:
                fault_matrix.loc[name,name2] = 1
    return fault_matrix