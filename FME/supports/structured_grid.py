import numpy as np
from .base_support import BaseGrid
class StructuredGrid(BaseGrid):
    def __init__(self,shape,step_vector=np.ones(3),origin=np.zeros(3),**kwargs):
        self.origin = origin
        self.shape = shape
        self.step_vector=step_vector
        self.n = shape[0]*shape[1]*shape[2]
        self.n_cell_x = self.shape[0]-1
        self.n_cell_y = self.shape[1] - 1
        self.n_cell_z = self.shape[2] - 1
        self.n_cell = self.n_cell_x*self.n_cell_y*self.n_cell_z
        self.aa = True
        if 'aa' in kwargs:
            self.aa = kwargs['aa']
        #TODO add in PCA for axis alignment
    def position_is_inside(self,pos):
        #check if point is inside mesh

    def position_to_dof_coefs(self,x,y):

