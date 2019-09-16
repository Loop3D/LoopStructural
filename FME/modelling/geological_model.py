from FME.modelling.geological_points import GPoint, IPoint
import numpy as np
import networkx as nx

class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes the
    """
    def __init__(self, origin, maximum, step_vector):
        """

        Parameters
        ----------
        origin - numpy array specifying the origin of the model
        maximum - numpy array specifying the maximum extent of the model
        """
        self.model = nx.Graph()
        self.features = {}
        self.data = {}
        self.data['gradient'] = []
        self.data['value'] = []
        self.origin = origin
        self.maximum = maximum
        self.step_vector = step_vector
        nsteps = (self.maximum-self.origin)/self.step_vector
        x = np.linspace(self.origin[0], nsteps[0] * self.step_vector[0], nsteps[0])
        y = np.linspace(self.origin[1], nsteps[1] * self.step_vector[1], nsteps[1])
        z = np.linspace(self.origin[2], nsteps[2] * self.step_vector[2], nsteps[2])
        self.xx, self.yy, self.zz = np.meshgrid(x, y, z, indexing='ij')

    def add_data(self, data):
        if type(data) == IPoint:
            self.data['value'].append(data)
        if type(data) == GPoint:
            self.data['gradient'].append(data)

    def add_feature(self, feature, name):
        self.features[name] = feature

    def evaluate(self):

        pass

    def plot_model_surface(self, featurename,**kwargs):

        pass

    def plot_model_scalar_field(self, featurename,**kwargs):

        pass

    def get_tetmesh_support(self,ntetra):

        pass

