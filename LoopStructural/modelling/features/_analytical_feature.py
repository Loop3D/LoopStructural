import numpy as np
from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
from typing import Optional

logger = getLogger(__name__)


class AnalyticalGeologicalFeature(BaseFeature):
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

    def __init__(self, name, vector, origin, regions=[], faults=[], model=None, builder=None):
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.vector = np.array(vector, dtype=float)
        self.origin = np.array(origin, dtype=float)
        self.type = FeatureType.ANALYTICAL

    def to_json(self):
        """
        Returns a json representation of the geological feature

        Returns
        -------
        json : dict
            json representation of the geological feature
        """
        json = super().to_json()
        json["vector"] = self.vector.tolist()
        json["origin"] = self.origin.tolist()
        return json

    def evaluate_value(self, xyz, ignore_regions=False):
        xyz = np.array(xyz)
        if len(xyz.shape) == 1:
            xyz = xyz[None, :]
        if len(xyz.shape) != 2:
            raise ValueError("xyz must be a 1D or 2D array")
        xyz2 = np.zeros(xyz.shape)
        xyz2[:] = xyz[:]
        for f in self.faults:
            xyz2[:] = f.apply_to_points(xyz)
        if self.model is not None:
            xyz2[:] = self.model.rescale(xyz2, inplace=False)
        xyz2[:] = xyz2 - self.origin
        normal = self.vector / np.linalg.norm(self.vector)
        distance = normal[0] * xyz2[:, 0] + normal[1] * xyz2[:, 1] + normal[2] * xyz2[:, 2]
        return distance / np.linalg.norm(self.vector)

    def evaluate_gradient(self, xyz, ignore_regions=False):
        xyz = np.array(xyz)
        if len(xyz.shape) == 1:
            xyz = xyz[None, :]
        if len(xyz.shape) != 2:
            raise ValueError("xyz must be a 1D or 2D array")
        v = np.zeros(xyz.shape)
        v[:, :] = self.vector[None, :]
        return v

    def get_data(self, value_map: Optional[dict] = None):
        return
