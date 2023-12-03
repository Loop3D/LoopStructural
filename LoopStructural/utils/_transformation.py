import numpy as np
from sklearn import decomposition
from sklearn.preprocessing import StandardScaler


class EuclideanTransformation:
    def __init__(self, dimensions=2):
        """Transforms points into a new coordinate
        system where the main eigenvector is aligned with x

        Parameters
        ----------
        dimensions : int, optional
            Do transformation in map view or on 3d volume, by default 2
        """
        self.rotation = None
        self.translation = None
        self.dimensions = dimensions
        self.angle = 0

    def fit(self, points: np.ndarray):
        """Fit the transformation to a point cloud
        This function will find the main eigenvector of the point cloud
        and rotate the point cloud so that this is aligned with x


        Parameters
        ----------
        points : np.ndarray
            xyz points as as numpy array
        """
        points = np.array(points)
        if points.shape[1] < self.dimensions:
            raise ValueError(
                "Points must have at least {} dimensions".format(self.dimensions)
            )
        # standardise the points so that centre is 0
        self.translation = np.mean(points, axis=0)
        # find main eigenvector and and calculate the angle of this with x
        pca = decomposition.PCA(n_components=self.dimensions).fit(
            points[:, : self.dimensions] - self.translation[None, : self.dimensions]
        )
        coeffs = pca.components_
        self.angle = -np.arccos(np.dot(coeffs[0, :], [1, 0]))
        self.rotation = self._rotation(self.angle)

    def _rotation(self, angle):
        return np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )

    def fit_transform(self, points: np.ndarray) -> np.ndarray:
        self.fit(points)
        return self.transform(points)

    def transform(self, points: np.ndarray) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        points : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        return np.dot(points - self.translation, self.rotation)

    def inverse_transform(self, points: np.ndarray) -> np.ndarray:
        return np.dot(points, self._rotation(-self.angle)) + self.translation

    def __call__(self, points: np.ndarray) -> np.ndarray:
        return self.transform(points)
