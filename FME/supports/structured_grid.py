import numpy as np
from .base_grid import BaseGrid
class StructuredGrid(BaseGrid):
    def __init__(self, shape=np.array([10,10,10]),step_vector=np.ones(3),origin=np.zeros(3),**kwargs):
        self.origin = origin
        self.shape = shape
        self.step_vector=step_vector
        self.n = shape[0]*shape[1]*shape[2]
        self.n_cell_x = self.shape[0]-1
        self.n_cell_y = self.shape[1] - 1
        self.n_cell_z = self.shape[2] - 1
        self.n_cell = self.n_cell_x*self.n_cell_y*self.n_cell_z
        self.coords = []
        for i in range(3):
            self.coords.append(np.linspace(self.origin[i],
                                  self.shape[i]*self.step_vector[0]+self.origin[0],
                                  self.shape[i]))

    def centers(self):
        xx, yy, zz = np.meshgrid(self.coords[0],
                                 self.coords[1],
                                 self.coords[2])
        xx = xx.flatten()
        yy = yy.flatten()
        zz = zz.flatten()
        return np.vstack([xx, yy, zz])
        #TODO add in PCA for axis alignment
    # def position_is_inside(self,pos):
    #     #check if point is inside mesh

    # def position_to_dof_coefs(self,x,y):

