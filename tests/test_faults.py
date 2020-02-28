from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_intrusion
data, bb = load_intrusion()
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)
fault = model.create_and_add_fault('fault',
                                   500,
                                   nelements=2000,
                                   steps=4,
                                   interpolatortype='PLI',
                                   buffer=0.2)
strati = model.create_and_add_foliation('strati',
                                        nelements=30000,
                                        interpolatortype='PLI',
                                        cgw=0.1,
                                        buffer=0.1
                                        )
