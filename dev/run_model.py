from LoopStructural.modelling.core.geological_model_graph import GeologicalModel
from LoopStructural.datasets import load_claudius

data, bb = load_claudius()

model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)
strati = model.create_and_add_foliation('strati')
strati.evaluate_value(model.regular_grid())