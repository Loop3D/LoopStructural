import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class AnalyticalGeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    model. For example foliations, fault planes, fold rotation angles etc. 
    
    Attributes
    ----------
    name : string
        should be a unique name for the geological feature
    support : a ScalarField
        holds the property values for the feature and links to the
        support geometry
    data : list
        list containing geological data
    region : list 
        list of boolean functions defining whether the feature is
        active
    faults : list 
        list of FaultSegments that affect this feature
    """
    def __init__(self, name, vector, origin, region=None, type=None, faults=[]):
        self.name = name
        self.vector = np.array(vector,dtype=float)
        self.origin = np.array(origin, dtype=float)
        self.region = region
        self.type = type
        self.faults = faults
        self.model = None

    def __str__(self):
        return self.name

    def __getitem__(self,key):
        return self._attributes[key]

    def __setitem__(self, key, item):
        self._attributes[key] = item

    def set_model(self, model):
        self.model = model

    def evaluate_value(self,xyz):
        xyz2 = np.zeros(xyz.shape)

        for f in self.faults:
            print('applying fault')
            xyz2[:] = f.apply_to_points(xyz)
        xyz2[:] = self.model.rescale(xyz2,inplace=False)
        xyz2[:] = xyz2-self.origin
        normal = self.vector/np.linalg.norm(self.vector)
        distance = normal[0]*xyz2[:,0]+normal[1]*xyz2[:,1]+normal[2]*xyz2[:,2]
        return distance/np.linalg.norm(self.vector)
    def evaluate_gradient(self,xyz):
        v = np.zeros(xyz.shape)
        v[:,:] = self.vector[None,:]
        return v 
    def min(self):
        if self.model is None:
                    return 0
        return np.nanmin(
                    self.evaluate_value(self.model.regular_grid((10, 10, 10))))
    def max(self):
        if self.model is None:
            return 0
        return np.nanmax(
            self.evaluate_value(self.model.regular_grid((10, 10, 10))))
