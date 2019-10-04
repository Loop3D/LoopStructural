import numpy as np

class BaseGrid:
    """
    A base class for a model support handles the coordinate transformations
    """
    def __init__(self, corners):
        """
        base constructor
        :param corners:
        """
        self.corners = corners
        

    def update_properties(self,propertyname,values):
        if values.shape[0] == self.n:
            self.proeprties[propertyname] = values
        if values.shape[0] == self.n_cell:
            self.cell_properties[propertyname] = values

    def is_inside(self, evaluation_points):
        transformed_evaluation_points = self.pca.transform(evaluation_points)
        inside = np.zeros(evaluation_points.shape[0]).astype(bool)
        inside[:] = True
        inside = transformed_evaluation_points[:, 0] > self.minpc0[None]
        inside *= transformed_evaluation_points[:, 0] < self.maxpc0[None]
        inside *= transformed_evaluation_points[:, 1] > self.minpc1[None]
        inside *= transformed_evaluation_points[:, 1] < self.maxpc1[None]
        inside *= transformed_evaluation_points[:, 2] > self.minpc2[None]
        inside *= transformed_evaluation_points[:, 2] < self.maxpc2[None]
        return inside

    def evaluate_value(self, evaluation_points, property_name):
        """
        virtual function for base class
        :param evaluation_points:
        :return:
        """
        pass

    def evaluate_gradient(self, evaluation_points, property_name):
        """
        virtual function for base class
        :param evaluation_points:
        :return:
        """

        pass

