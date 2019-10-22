from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
import pandas as pd
import numpy as np

model = GeologicalModel((0,0,0),(4,4,4))
df = pd.read_csv('model_input.csv')
model.set_model_data(df)
s0 = model.create_and_add_conformable_foliation('s0',
                                                interpolatortype='FDI',
                                                nelements=20000,
                                                solver='cg',
                                                damp=False)

# uc = model.create_and_add_unconformity('uc',
#                                        interpolatortype='FDI',
#                                        nelements=20000,
#                                        solver='cg',
#                                        damp=False)
# fault = model.create_and_add_fault('fault', .3,
#                                       interpolatortype='FDI',
#                                       nelements=20000,
#                                       solver='cg',
#                                       damp=False)
# #
s0_2 = model.create_and_add_conformable_foliation('s0_2',
                                                  interpolatortype='FDI',
                                                  nelements=20000,
                                                  solver='cg',
                                                  damp=False)
viewer = LavaVuModelViewer(background='white')
try:
    viewer.add_isosurface(s0_2,
                      # slices=[0,1,2],#3#isovalue=1,
                      nslices=10,
                    voxet=model.voxet(),
                      colour='black'
                  )
except:
    pass
try:
    viewer.add_isosurface(fault.faultframe[0],
                          isovalue=0,
                          # nslices=4,
                        voxet=model.voxet(),
                          colour='black'
                      )
except:
    print('fault')
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
# viewer.add_data(uc)
# viewer.add_scalar_field(model.bounding_box,(50,50,50),
#                           'scalar',
# #                             norm=True,
#                          paint_with=s0_2,
#                          cmap='tab20')
# print(uc.get_interpolator().get_value_constraints())
# print(uc.get_node_values())
# print(uc.evaluate_value(uc.get_interpolator().get_value_constraints()[:,:3]))
viewer.interactive()
