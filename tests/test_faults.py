from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
# from LoopStructural.datasets import load_intrusion
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)
fault = model.create_and_add_fault('fault',
                                   500,
                                   nelements=2000,
                                   steps=4,
                                   interpolatortype='PLI')
strati = model.create_and_add_foliation('strati',nelements=30000,interpolatortype='PLI',cgw=0.1)