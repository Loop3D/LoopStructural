from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
from LoopStructural import GeologicalModel
from LoopStructural.utils import log_to_file
import numpy as np
# log_to_file('dev/ls.log')


fault_params = {'interpolatortype':'FDI',
                'nelements':1e3,
                'data_region':.1000, 
                'solver':'pyamg',
#                 overprints:overprints,
                'cpw':10,
                'npw':10}
foliation_params = {'nelements':1e4,  # how many tetras/voxels
                    'buffer':2.5,  # how much to extend nterpolation around box
                    'solver':'pyamg',
                    'npw':10,
                    'cpw':10,
                    'interpolatortype':'FDI',
                    'damp':True}

model, m2l_data = GeologicalModel.from_map2loop_directory('./dev/unconf',
                                                        skip_faults=False,
                                                          rescale=False,
                                                        fault_params=fault_params,
                                                        foliation_params=foliation_params)
view = LavaVuModelViewer(model,vertical_exaggeration=1) 
for f in model['supergroup_0'].faults:
  bb = np.zeros((3,2))
  bb[:,0] = f.builder.origin
  bb[:,1] = f.builder.maximum
  view.add_box(bb,f.name+'_support')
view.interactive()
# print(model.data[model.data['coord']==1])
# # print(model.features[1][0].builder.data)
# for f in features:
#     f.builder.update()
# model.features[1][0].builder.update()#evaluate_value(np.array([[0,0,0]]))
# print(model.features[1][0].interpolator.support)
# print(model.support)
# print(model['supergroup_0'].faults[0][0].builder.build_arguments)
# model['supergroup_0'].faults[0][0].builder.update()
# print(model)
# for f in model['supergroup_0'].faults:
#     print(f[0].interpolator.nx)
# view = LavaVuModelViewer(model,vertical_exaggeration=1) 
# # for i in range(3):
# #     view.add_isosurface(model['supergroup_0'].faults[i][0],value=0)
# # view.interactive()
# # for f in model['supergroup_0'].faults:'
# #     # f = model['supergroup_0'].faults[i]
# #     view.add_vector_field(f,locations=model.regular_grid()[::100,:])
# # # # model['supergroup_0'].evaluate_value(model.regular_grid(nsteps=(10,10,10)))
# #     view.add_isosurface(f,value=0)
# # # # view.interactive()
# view.add_model_surfaces()#filename=filename)
# # # # # view.add_data(model['supergroup_0'])
# # # view.add_fault_displacements()
# view.interactive()