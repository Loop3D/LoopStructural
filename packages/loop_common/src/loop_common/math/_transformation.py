from __future__ import annotations

import numpy as np

from loop_common.logging import get_logger

logger = get_logger(__name__)


def _import_decomposition():
    """Import the scikit-learn decomposition module.

    Returns
    -------
    module
        The module ``sklearn.decomposition``.

    Raises
    ------
    ImportError
        If scikit-learn is not installed.
    """
    try:
        from sklearn import decomposition
    except ImportError as error:
        raise ImportError(
            "scikit-learn is required for this function. Install it with: pip install scikit-learn"
        ) from error
    return decomposition


class EuclideanTransformation:
    def __init__(
        self,
        dimensions: int = 2,
        angle: float = 0,
        translation: np.ndarray | None = None,
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
        if translation is None:
            translation = np.zeros(3)
        self.translation = translation[:dimensions]
        self.dimensions = dimensions
        self.angle = angle
        self.fit_rotation = fit_rotation

    def _check_points(self, points) -> np.ndarray:
        """Check the points and return a float copy.

        Parameters
        ----------
        points : array_like
            Points with shape (N, d).

        Returns
        -------
        np.ndarray
            A copy of the points as a float array.

        Raises
        ------
        ValueError
            If the array is not 2D, or has fewer columns than ``dimensions``.
        """
        array = np.array(points, dtype=float)
        if array.ndim != 2:
            raise ValueError("points must be a 2D array with shape (N, d)")
        if array.shape[1] < self.dimensions:
            raise ValueError(f"Points must have at least {self.dimensions} dimensions")
        return array

    def fit(self, points: np.ndarray):
        """Fit the transformation to a point cloud.

        Find the centre of the points. If ``fit_rotation`` is true, find the
        main axis of the points in the xy plane and set the angle so that
        this axis is aligned with x. With fewer than 2 rows, the angle is 0
        and scikit-learn is not used. The object does not change if this
        method raises an error.

        Parameters
        ----------
        points : np.ndarray
            Points with shape (N, d). d must be at least ``dimensions``.

        Returns
        -------
        EuclideanTransformation
            The fitted object (self).

        Raises
        ------
        ValueError
            If the points are not valid.
        ImportError
            If ``fit_rotation`` is true, N is 2 or more, and scikit-learn
            is not installed.
        """
        array = self._check_points(points)
        d = self.dimensions
        translation = np.mean(array[:, :d], axis=0)
        angle = 0
        # PCA needs at least 2 rows. With 1 row the angle stays 0.
        if self.fit_rotation and array.shape[0] >= 2:
            decomposition = _import_decomposition()
            centred = array[:, :2] - translation[None, :2]
            axis = decomposition.PCA(n_components=2).fit(centred).components_[0]
            angle = -np.arctan2(axis[1], axis[0])
        self.translation = translation
        self.angle = angle
        logger.debug(f"Fitted transformation: translation {translation}, angle {angle}")
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
                [0, 0, 1],
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
        array = self._check_points(points)
        d = self.dimensions
        centred = array[:, :d] - self.translation[None, :]
        array[:, :d] = np.einsum("ik,jk->ij", centred, self.rotation[:d, :d])
        return array

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
        array = self._check_points(points)
        d = self.dimensions
        array[:, :d] = (
            np.einsum("ik,jk->ij", array[:, :d], self.inverse_rotation[:d, :d])
            + self.translation[None, :]
        )
        return array

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
        html_str = f"""
        <div class="collapsible">
          <button class="collapsible-button">{self.__class__.__name__}</button>
          <div class="content">
            <p>Translation: {self.translation}</p>
            <p>Rotation Angle: {self.angle} degrees</p>
          </div>
        </div>
        """
        return html_str
