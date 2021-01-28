# from LoopStructural.modelling.core.geological_model_graph import GeologicalModel
from LoopStructural import GeologicalModel
import numpy as np
from LoopStructural.datasets import load_intrusion, load_claudius
from LoopStructural.utils import log_to_file
from LoopStructural.visualisation import LavaVuModelViewer

data, bb = load_intrusion()
log_to_file('test.log')
model = GeologicalModel(bb[0,:],bb[1,:],rescale=False)

model.set_model_data(data)
fault = model.create_and_add_fault('fault',
                                    600,
                                    interpolatortype='FDI',
                                    nelements=1e3,
                                    solver='pyamg',
                                    steps=2
                                    )
strati = model.create_and_add_foliation('strati',
                                        interpolatortype='FDI',
                                        solver='pyamg')
print(model )
strati.evaluate_value(model.regular_grid())
# # print(strati.faults)
# print(fault[0].faults)
# view = LavaVuModelViewer(model,background='white')
# view.add_isosurface(strati,slices=np.unique(data.loc[~np.isnan(data['val']),'val']))
# view.add_isosurface(fault,value=0)
# view.add_data(strati)
# # view.add_vector_field(fault,locations=model.regular_grid()[::10,:])
# view.interactive()