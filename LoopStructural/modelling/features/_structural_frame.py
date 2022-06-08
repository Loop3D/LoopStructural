"""
Structural frames
"""
import logging
from LoopStructural.modelling.features import BaseFeature
import numpy as np
from LoopStructural.utils import getLogger
from LoopStructural.modelling.features import FeatureType

logger = getLogger(__name__)


class StructuralFrame(BaseFeature):
    def __init__(self, name, features, fold=None):
        """
        Structural frame is a curvilinear coordinate system defined by
        structural observations associated with a fault or fold.

        Parameters
        ----------
        name - name of the structural frame
        features - list of features to build the frame with
        """
        BaseFeature.__init__(self, name, None, [], [], None)
        self.features = features
        self.fold = fold
        self.type = FeatureType.STRUCTURALFRAME

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
        for i in range(3):
            self.features[i].add_region(region)

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

    def evaluate_value(self, evaluation_points):
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
        v = np.zeros(evaluation_points.shape)  # create new 3d array of correct length
        v[:] = np.nan
        v[:, 0] = self.features[0].evaluate_value(evaluation_points)
        v[:, 1] = self.features[1].evaluate_value(evaluation_points)
        v[:, 2] = self.features[2].evaluate_value(evaluation_points)
        return v

    def evaluate_gradient(self, evaluation_points, i=None):
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
            return self.features[i].support.evaluate_gradient(evaluation_points)
        return (
            self.features[0].support.evaluate_gradient(evaluation_points),
            self.features[1].support.evaluate_gradient(evaluation_points),
            self.features[2].support.evaluate_gradient(evaluation_points),
        )
