import logging

import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class FoldRotationAngleFeature:
    """

    """
    def __init__(self, fold_frame, rotation):
        """

        Parameters
        ----------
        fold_frame
        rotation
        """
        self.fold_frame = fold_frame
        self.rotation = rotation

    def evaluate_value(self, location):
        """

        Parameters
        ----------
        location

        Returns
        -------

        """
        s1 = self.fold_frame.features[0].evaluate_value(location)
        r = self.rotation(s1)
        return r


def fourier_series(x, c0, c1, c2, w):
    """

    Parameters
    ----------
    x
    c0
    c1
    c2
    w

    Returns
    -------

    """
    v = np.array(x.astype(float))
    # v.fill(c0)
    v = c0 + c1 * np.cos(2 * np.pi / w * x) + c2 * np.sin(2 * np.pi / w * x)
    return np.rad2deg(np.arctan(v))
