from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
import pandas as pd
import numpy as np

model = GeologicalModel((-4,-4,-4),(4,4,4))
print(model.bounding_box)
df = pd.read_csv('model_input.csv')
model.set_model_data(df)
s0_2 = model.create_and_add_conformable_series('s0_2', solver='cg', damp=False)
uc = model.create_and_add_unconformity('uc', solver='cg', damp=False)
s0 = model.create_and_add_conformable_series('s0', solver='cg', damp=False)

viewer = LavaVuModelViewer(background='white')
try:
    viewer.add_isosurface(s0_2,
                          isovalue=1,
                          # nslices=4,
                        voxet=model.voxet(),
                          colour='black'
                      )
except:
    print('so2')
    pass
try:
    viewer.add_isosurface(s0,
                        voxet=model.voxet(),
                          nslices=10,
                          colour='green'

                          )
except:
    print('s0')
    pass
try:
    viewer.add_isosurface(uc,
                      slices=[0],
                        voxet=model.voxet()
                      )
except:
    print('uc')
    pass
viewer.add_data(uc)
# viewer.add_scalar_field(model.bounding_box,(50,10,50),
#                           'scalar',
# #                             norm=True,
#                          paint_with=uc,
#                          cmap='tab20')
print(uc.get_interpolator().get_value_constraints())
print(uc.get_node_values())
print(uc.evaluate_value(uc.get_interpolator().get_value_constraints()[:,:3]))
viewer.interactive()