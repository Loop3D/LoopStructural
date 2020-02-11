import numpy as np
from LoopStructural.supports.structured_tetra import TetMesh
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.utils.helper import strike_dip_vector, plunge_and_plunge_dir_to_vector

mesh = TetMesh(nsteps=[3,3,3])

# print('els',mesh.get_elements())
neighbours2 = mesh.get_neighbours()
# print('mask1',mesh.tetra_mask)
# print('mask2',mesh.tetra_mask_even)
# print(mesh.n_cell_x,mesh.n_cell_y,mesh.n_cell_z)
neighbours = np.zeros((mesh.n_elements,4)).astype(int)
neighbours[:] = -1
elements = mesh.get_elements()
print('ee',elements)
for ie, e in enumerate(elements):
    nn = 0
    for iin, ne in enumerate(elements):
        n = 0
        for i in range(4):
            for j in range(4):
                if e[i] == ne[j]:
                    n+=1
        if n == 3:
            neighbours[ie,nn] = iin
            nn+=1
for i in range(neighbours.shape[0]):
    for j in range(4):
        if neighbours[i,j] not in neighbours2[i,:]:
            print(neighbours[i,:],neighbours2[i,:])
    #   print(neighbours[i,:],neighbours2[i,:])
