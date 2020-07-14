"""
Structural frames
"""
import logging

logger = logging.getLogger(__name__)


class StructuralFrame:
    """[summary]

    [extended_summary]
    """
    def __init__(self, name, features, fold=None):
        """
        Structural frame is a curvilinear coordinate system defined by
        structural
        observations associated with a fault or fold.

        Parameters
        ----------
        name - name of the structural frame
        features - list of features to build the frame with
        """
        self.name = name
        self.features = features
        self.data = None
        self.fold = fold
    def __getitem__(self, item):
        """

        Parameters
        ----------
        item index of feature to access

        Returns
        -------
        the structural frame geological feature
        """
        return self.features[item]

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

    def set_data(self, data):
        """
        Associate data with structural frame

        Parameters
        ----------
        data

        Returns
        -------

        """
        self.data = data

    def evaluate_value(self, evaluation_points, i=None):
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
        if i is not None:
            self.features[i].support.evaluate_value(evaluation_points)
        return (self.features[0].support.evaluate_value(evaluation_points),
                self.features[1].support.evaluate_value(evaluation_points),
                self.features[2].support.evaluate_value(evaluation_points))

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
            return self.features[i].support.evaluate_gradient(
                evaluation_points)
        return (self.features[0].support.evaluate_gradient(evaluation_points),
                self.features[1].support.evaluate_gradient(evaluation_points),
                self.features[2].support.evaluate_gradient(evaluation_points))



