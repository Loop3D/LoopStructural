from LoopStructural.modelling.core.geological_model_graph import GeologicalModel
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
                                    solver='pyamg',
                                    steps=2
                                    )
strati = model.create_and_add_foliation('strati',
                                        interpolatortype='FDI',
                                        solver='pyamg')

# strati.evaluate_value(model.regular_grid())
print(strati.faults)
print(fault[0].faults)
view = LavaVuModelViewer(model)
view.add_isosurface(strati)
view.add_isosurface(fault)
view.interactive()