import logging

import numpy as np

logger = logging.getLogger(__name__)


class ScalarField:
    def __init__(self, support, property_name):
        """
         A scalar field is a distance from a reference location/horizon
         within a model area
        currently the scalar field just interfaces between the model support
        and the geological
        feature. In the future, the property values might be stored on the
        scalar field rather
        than as a dict on the tetmesh?
        Parameters
        ----------
        support -
            some geometrical object that represents a discrete supprort such
            as a
            tetmesh or structured grid
        property_name - string
            name of the property saved on the support
        """
        self.support = support
        self.property_name = property_name
        self.interpolator = None

    @classmethod
    def from_node_values(cls, support, property_name, node_values):
        """
        Build a scalar field from an array of node values on a support
        saving the
        values in the property dictionary
        Parameters
        ----------
        support
        property_name
        node_values

        Returns
        -------

        """
        support.update_property(property_name, node_values)
        return cls(support, property_name)

    @classmethod
    def from_interpolator(cls, interpolator):
        """
        Buidl a scalar field using an interpolator
        Parameters
        ----------
        interpolator

        Returns
        -------

        """
        interpolator.update()
        interpolator.support.update_property(interpolator.propertyname,
                                             interpolator.c)
        scalar_field = cls(
            interpolator.support,
            interpolator.propertyname)
        scalar_field.interpolator = interpolator
        return scalar_field

    def evaluate_value(self, evaluation_points):
        """
        Evaluate the scalar value fo the scalar field at the locations
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_value(
                evaluation_points[~mask], self.property_name)
        return evaluated

    def evaluate_gradient(self, evaluation_points):
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        if evaluation_points.shape[0] > 0:
            return self.support.evaluate_gradient(evaluation_points,
                                                  self.property_name)
        return np.zeros((0, 3))

    def mean(self):
        """
        average value of the scalar field
        Returns
        -------

        """
        return np.nanmean(self.support.properties[self.property_name])

    def min(self):
        """
        Min value of the scalar field
        Returns
        -------

        """
        return np.nanmin(self.support.properties[self.property_name])

    def max(self):
        return np.nanmax(self.support.properties[self.property_name])

    def get_node_values(self):
        """
        Node values from the mesh object
        Returns
        -------

        """
        return self.support.properties[self.property_name]

    def number_of_nodes(self):
        """
        Number of nodes in the mesh
        Returns
        -------

        """
        return self.support.n

    def update_property(self, values):
        """
        Updates the values of the scalar field on the mesh
        Parameters
        ----------
        values

        Returns
        -------

        """
        self.support.properties[self.property_name] = values

    def slice(self, isovalue, region):
        """
        Extract an isosurface from a scalar field
        Parameters
        ----------
        isovalue
        region

        Returns
        -------

        """
        return self.support.slice(self.property_name, isovalue, "everywhere")
