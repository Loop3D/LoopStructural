from LoopStructural.interpolators import StructuredGrid
import numpy as np
class ScalarField:
    """A scalar field defined by a regular grid and values
    """
    def __init__(self,values,step_vector,nsteps, name='scalar field'):
        """[summary]

        Parameters
        ----------
        values : numpy array
            nodes values of regular grid
        step_vector : numpy array
            step vector for x,y,z step
        nsteps : numpy array
            number of steps
        name : string, optional
            name of the feature for the visualisation
        """
        self.values = values
        self.grid = StructuredGrid(np.zeros(3),nsteps,step_vector)
        self.name = name
        
    @property
    def nodes(self):
        return self.grid.nodes


    def evaluate_value(self,xyz):
        """Evaluate the scalar field at locations

        Parameters
        ----------
        xyz : numpy array
            locations in real coordinates

        Returns
        -------
        numpy array
            interpolated values
        """
        v = self.grid.evaluate_value(xyz,self.values)
        return v

    def min(self):
        return np.min(self.values)
    def max(self):
        return np.max(self.values)#14