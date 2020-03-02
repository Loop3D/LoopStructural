from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
from LoopStructural.datasets import load_intrusion
import pandas as pd
import numpy as np
data = pd.read_pickle('../../notebooks/scratch/newdata.pkl')
bb = np.loadtxt('../../notebooks/scratch/bb.txt')
# print(bb)
bb[0,2]=-5000
bb[1,2] = 5000
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)

# fault = model.create_and_add_fault('fault',-500,nelements=5000,steps=4,interpolatortype='PLI',buffer=0.2)
strati = model.create_and_add_foliation('fault',
                                        nelements=100,
                                        interpolatortype='PLI',
                                        regularisation=0.3,
                                        buffer=0.5
                                       )

# mesh = strati['feature'].get_interpolator().support
# print(mesh.origin,mesh.nsteps_cells,mesh.step_vector)
# mesh.max = mesh.origin + mesh.nsteps_cells * mesh.step_vector
# x = np.linspace(mesh.origin[0], mesh.max[0], mesh.nsteps[0])
# y = np.linspace(mesh.origin[1], mesh.max[1], mesh.nsteps[1])
# z = np.linspace(mesh.origin[2], mesh.max[2], mesh.nsteps[2])
# xx, yy, zz = np.meshgrid(x, y, z,indexing='ij')
# mesh.nodes = np.array([xx.flatten(order='F'), yy.flatten(order='F'),
#                        zz.flatten(order='F')]).T
#
# strati['feature'].get_interpolator().support.properties['fault'] = \
#                 strati['feature'].get_interpolator().support.nodes[:,2]
viewer = LavaVuModelViewer(model)
viewer.add_scalar_field(strati['feature'],locations=model.regular_grid()[::40])
# viewer.add_isosurface(strati['feature'])
# viewer.add_data(strati['feature'])
# viewer.add_points(model.data[model.data['type']=='strati'][['X','Y','Z']],name='prefault',pointsize=5,colour='blue')
viewer.interactive()