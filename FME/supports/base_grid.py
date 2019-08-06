from sklearn.decomposition import PCA

3
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
        self.pca = PCA(n_components=3)
        self.pca.fit(corners)
        self.transformed_corners = self.pca.transform(corners)
        self.minpc0 = np.min(self.transformed_corners[:, 0])
        self.maxpc0 = np.max(self.transformed_corners[:, 0])
        self.minpc1 = np.min(self.transformed_corners[:, 1])
        self.maxpc1 = np.max(self.transformed_corners[:, 1])
        self.minpc2 = np.min(self.transformed_corners[:, 2])
        self.maxpc2 = np.max(self.transformed_corners[:, 2])

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

