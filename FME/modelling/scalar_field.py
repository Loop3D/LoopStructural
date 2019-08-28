import numpy as np


class ScalarField:
    def __init__(self, support, property_name):
        self.support = support
        self.property_name = property_name
        self.interpolator = None

    @classmethod
    def from_node_values(cls, support, property_name, node_values):
        support.update_property(property_name, node_values)
        return cls(support, property_name)

    @classmethod
    def from_interpolator(cls, interpolator):
        interpolator.update()
        interpolator.support.update_property(interpolator.propertyname, interpolator.c)
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
        mask = np.any(np.isnan(evaluation_points),axis=1)

        if evaluation_points[~mask,:].shape[0]>0:
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
        if evaluation_points.shape[0]>0:
            return self.support.evaluate_gradient(evaluation_points, self.property_name)
        return np.zeros((0,3))

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

    def slice(self, isovalue):
        return self.support.slice(self.property_name, isovalue)