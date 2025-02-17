import numpy as np
from . import getLogger

logger = getLogger(__name__)


class EuclideanTransformation:
    def __init__(
        self,
        dimensions: int = 2,
        angle: float = 0,
        translation: np.ndarray = np.zeros(3),
        fit_rotation: bool = True,
    ):
        """Transforms points into a new coordinate
        system where the main eigenvector is aligned with x

        Parameters
        ----------
        dimensions : int, optional
            Do transformation in map view or on 3d volume, by default 2
        angle : float, optional
            Angle to rotate the points by, by default 0
        translation : np.ndarray, default zeros
            Translation to apply to the points, by default
        """
        self.translation = translation[:dimensions]
        self.dimensions = dimensions
        self.angle = angle
        self.fit_rotation = fit_rotation

    def fit(self, points: np.ndarray):
        """Fit the transformation to a point cloud
        This function will find the main eigenvector of the point cloud
        and rotate the point cloud so that this is aligned with x


        Parameters
        ----------
        points : np.ndarray
            xyz points as as numpy array
        """
        try:
            from sklearn import decomposition
        except ImportError:
            logger.error('scikit-learn is required for this function')
            return
        points = np.array(points)
        if points.shape[1] < self.dimensions:
            raise ValueError("Points must have at least {} dimensions".format(self.dimensions))
        # standardise the points so that centre is 0
        # self.translation = np.zeros(3)
        self.translation = np.mean(points[:, : self.dimensions], axis=0)
        # find main eigenvector and and calculate the angle of this with x
        if self.fit_rotation:
            pca = decomposition.PCA(n_components=self.dimensions).fit(
                points[:, : self.dimensions] - self.translation[None, : self.dimensions]
            )
            coeffs = pca.components_
            self.angle = -np.arccos(np.dot(coeffs[0, :], [1, 0]))
        else:
            self.angle = 0
        return self

    @property
    def rotation(self):
        return self._rotation(self.angle)

    @property
    def inverse_rotation(self):
        return self._rotation(-self.angle)

    def _rotation(self, angle):
        return np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, -1],
            ]
        )

    def fit_transform(self, points: np.ndarray) -> np.ndarray:
        """Fit the transformation and transform the points"""

        self.fit(points)
        return self.transform(points)

    def transform(self, points: np.ndarray) -> np.ndarray:
        """Transform points using the transformation and rotation

        Parameters
        ----------
        points : np.ndarray
            xyz points as as numpy array

        Returns
        -------
        np.ndarray
            xyz points in the transformed coordinate system
        """
        points = np.array(points)
        if points.shape[1] < self.dimensions:
            raise ValueError("Points must have at least {} dimensions".format(self.dimensions))
        centred = points[:, : self.dimensions] - self.translation[None, :]
        rotated = np.einsum(
            'ik,jk->ij',
            centred,
            self.rotation[: self.dimensions, : self.dimensions],
        )
        transformed_points = np.copy(points)
        transformed_points[:, : self.dimensions] = rotated
        return transformed_points

    def inverse_transform(self, points: np.ndarray) -> np.ndarray:
        """
        Transform points back to the original coordinate system

        Parameters
        ----------
        points : np.ndarray
            xyz points as as numpy array

        Returns
        -------
        np.ndarray
            xyz points in the original coordinate system
        """
        inversed = (
            np.einsum(
                'ik,jk->ij',
                points[: self.dimensions],
                self.inverse_rotation[: self.dimensions, : self.dimensions],
            )
            + self.translation
        )
        inversed = (
            np.vstack([inversed, points[self.dimensions :]])
            if points.shape[1] > self.dimensions
            else inversed
        )
        return inversed

    def __call__(self, points: np.ndarray) -> np.ndarray:
        """
        Transform points into the transformed space

        Parameters
        ----------
        points : np.ndarray
            xyz points as as numpy array

        Returns
        -------
        np.ndarray
            xyz points in the transformed coordinate system
        """

        return self.transform(points)

    def _repr_html_(self):
        """
        Provides an HTML representation of the TransRotator.
        """
        html_str = """
        <div class="collapsible">
          <button class="collapsible-button">{self.__class__.__name__}</button>
          <div class="content">
            <p>Translation: {self.translation}</p>
            <p>Rotation Angle: {self.angle} degrees</p>
          </div>
        </div>
        """.format(
            self=self
        )
        return html_str
