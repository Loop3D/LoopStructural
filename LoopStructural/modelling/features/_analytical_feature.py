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

    def __init__(
        self,
        name: str,
        vector: np.ndarray,
        origin: np.ndarray,
        regions=[],
        faults=[],
        model=None,
        builder=None,
    ):
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        try:
            self.vector = np.array(vector, dtype=float).reshape(3)
        except ValueError:
            logger.error("AnalyticalGeologicalFeature: vector must be a 3 element array")
            raise ValueError("vector must be a 3 element array")
        try:
            self.origin = np.array(origin, dtype=float).reshape(3)
        except ValueError:
            logger.error("AnalyticalGeologicalFeature: origin must be a 3 element array")
            raise ValueError("origin must be a 3 element array")
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

    def evaluate_value(self, pos: np.ndarray, ignore_regions=False):
        pos = np.array(pos)
        if len(pos.shape) == 1:
            pos = pos[None, :]
        if len(pos.shape) != 2:
            raise ValueError("xyz must be a 1D or 2D array")
        xyz2 = np.zeros(pos.shape)
        xyz2[:] = pos[:]
        for f in self.faults:
            xyz2[:] = f.apply_to_points(pos)
        if self.model is not None:
            xyz2[:] = self.model.rescale(xyz2, inplace=False)
        xyz2[:] = xyz2 - self.origin
        normal = self.vector / np.linalg.norm(self.vector)
        distance = normal[0] * xyz2[:, 0] + normal[1] * xyz2[:, 1] + normal[2] * xyz2[:, 2]
        return distance / np.linalg.norm(self.vector)

    def evaluate_gradient(self, pos: np.ndarray, ignore_regions=False):
        pos = np.array(pos)
        if len(pos.shape) == 1:
            pos = pos[None, :]
        if len(pos.shape) != 2:
            raise ValueError("pos must be a 1D or 2D array")
        v = np.zeros(pos.shape)
        v[:, :] = self.vector[None, :]
        return v

    def get_data(self, value_map: Optional[dict] = None):
        return

    def copy(self, name: Optional[str] = None):
        if name is None:
            name = self.name
        return AnalyticalGeologicalFeature(
            name,
            self.vector.copy(),
            self.origin.copy(),
            list(self.regions),
            list(self.faults),
            self.model,
            self.builder,
        )
