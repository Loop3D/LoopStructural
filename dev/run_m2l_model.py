from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
from LoopStructural import GeologicalModel
from LoopStructural.utils import log_to_file
import numpy as np
# log_to_file('dev/ls.log')
fault_params = {'interpolatortype':'FDI',
                'nelements':1e4,
                # 'data_region':0,
                'solver':'pyamg',
                'step':10,
                'cpw':10,
                'npw':10,
                'tol':1e-8
               }
foliation_params = {'interpolatortype':'FDI' , # 'interpolatortype':'PLI',
                'nelements':1e4,  # how many tetras/voxels
                'data_region':.3,  # how much to extend nterpolation around box
                'solver':'pyamg',
                'damp':True}

fault_params = {'interpolatortype':'FDI',
                'nelements':1e4,
                # 'data_region':.2, 
                'solver':'lu',
#                 overprints:overprints,
                'cpw':10,
                'npw':10}
foliation_params = {'nelements':1e5,  # how many tetras/voxels
                    'buffer':2.5,  # how much to extend nterpolation around box
                    'solver':'pyamg',
                    'npw':10,
                    'cpw':10,
                    'interpolatortype':'FDI',
                    'damp':True}

model, m2l_data = GeologicalModel.from_map2loop_directory('./dev/unconf',
                                                        skip_faults=False,
                                                        #   rescale=False,
                                                        fault_params=fault_params,
                                                        foliation_params=foliation_params)
# print(model.data[model.data['coord']==1])
# # print(model.features[1][0].builder.data)
# for f in features:
#     f.builder.update()
# model.features[1][0].builder.update()#evaluate_value(np.array([[0,0,0]]))
# print(model.features[1][0].interpolator.support)
# print(model.support)
view = LavaVuModelViewer(model,vertical_exaggeration=1) 
view.add_isosurface(model.features[1][0])
# # view.interactive()
view.add_model_surfaces(faults=True)#filename=filename)
# # # view.add_data(model['supergroup_0'])
# view.add_fault_displacements()
view.interactive()