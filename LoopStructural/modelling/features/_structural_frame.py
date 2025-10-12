"""
Structural frames
"""

from ..features import BaseFeature, FeatureType
import numpy as np
from ...utils import getLogger
from typing import Optional, List, Union
from ...datatypes import ValuePoints, VectorPoints

logger = getLogger(__name__)


class StructuralFrame(BaseFeature):
    def __init__(self, name: str, features: list, fold=None, model=None):
        """
        Structural frame is a curvilinear coordinate system defined by
        structural observations associated with a fault or fold.

        Parameters
        ----------
        name - name of the structural frame
        features - list of features to build the frame with
        """
        BaseFeature.__init__(self, name, model, [], [], None)
        self.features = features
        self.fold = fold
        self.type = FeatureType.STRUCTURALFRAME

    def to_json(self):
        """
        Return a json representation of the structural frame

        Returns
        -------
        json : dict
            json representation of the structural frame
        """
        json = {}
        json["name"] = self.name
        json["features"] = [f.name for f in self.features]
        json["type"] = self.type
        return json

    def __getitem__(self, key):
        """

        Parameters
        ----------
        key index of feature to access

        Returns
        -------
        the structural frame geological feature
        """
        return self.features[key]

    def __setitem__(self, key, value):
        """

        Parameters
        ----------
        item index of feature to access

        Returns
        -------
        the structural frame geological feature
        """
        self.features[key] = value

    def set_model(self, model):
        """Link the model that created the frame to the frame
        and the features that make up the frame

        Parameters
        ----------
        model : GeologicalModel
            the geological model that created the fold frame
        """
        self.model = model
        for f in self.features:
            if f is None:
                continue
            f.set_model(model)

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        # causes circular import, could delay import?
        # if type(model) == GeologicalModel:
        for f in self.features:
            if f is None:
                continue
            f.model = model

    def add_region(self, region):
        self.regions.append(region)
        for i in range(3):
            self.features[i].regions = self.regions

    def get_feature(self, i):
        """
        Return the ith feature

        Parameters
        ----------
        i

        Returns
        -------

        """
        return self.features[i]

    def evaluate_value(self, pos, ignore_regions=False):
        """
        Evaluate the value of the structural frame for the points.
        Can optionally only evaluate one coordinate

        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        v = np.zeros(pos.shape)  # create new 3d array of correct length
        v[:] = np.nan
        v[:, 0] = self.features[0].evaluate_value(pos, ignore_regions=ignore_regions)
        v[:, 1] = self.features[1].evaluate_value(pos, ignore_regions=ignore_regions)
        v[:, 2] = self.features[2].evaluate_value(pos, ignore_regions=ignore_regions)
        return v[:,0]

    def evaluate_gradient(self, pos, i=None, ignore_regions=False):
        """
        Evaluate the gradient of the structural frame.
        Can optionally only evaluate the ith coordinate

        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        if i is not None:
            return self.features[i].interpolator.evaluate_gradient(pos)
        return self.features[0].interpolator.evaluate_gradient(pos)            

    def get_data(self, value_map: Optional[dict] = None) -> List[Union[ValuePoints, VectorPoints]]:
        """Return the data associated with the features in the
        structural frame

        Parameters
        ----------
        value_map : Optional[dict], optional
            map scalar values to another property, by default None

        Returns
        -------
        List
            container of value or vector points
        """
        data = []
        for f in self.features:
            data.extend(f.get_data(value_map))
        return data

    def copy(self, name: Optional[str] = None):
        if name is None:
            name = f'{self.name}_copy'
        # !TODO check if this needs to be a deep copy
        return StructuralFrame(name, self.features, self.fold, self.model)
