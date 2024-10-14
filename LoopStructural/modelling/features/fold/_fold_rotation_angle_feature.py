import numpy as np
from ....modelling.features import BaseFeature
<<<<<<< HEAD
from abc import ABCMeta, abstractmethod
=======
>>>>>>> master
from ....utils import getLogger

logger = getLogger(__name__)


class FoldRotationAngleFeature(BaseFeature):
    """ """

    def __init__(
        self,
        fold_frame,
        rotation,
        name="fold_rotation_angle",
        model=None,
        faults=[],
        regions=[],
        builder=None,
    ):
        """

        Parameters
        ----------
        fold_frame
        rotation
        """
        BaseFeature.__init__(self, f"{name}_displacement", model, faults, regions, builder)
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


# class BaseFoldProfile(ABCMeta):
#     def __init__(self):
#         pass

#     @abstractmethod
#     def params(self):
#         pass


#     @abstractmethod
#     def __call__(self, s):
#         pass


# class TrigoFoldProfile(BaseFoldProfile):
#     def __init__(self, origin, wavelength, inflectionpointangle):
#         self.origin = origin
#         self.wavelength = wavelength
#         self.inflectionpointangle = inflectionpointangle
#         self.tan_alpha_delta_half = np.tan(self.inflectionpointangle)

#         self.tan_alpha_shift = 0

#     def position_in_period(self, s):
#         return (s - self.origin) / self.wavelength

#     def __call__(self, s):
#         x = self.position_in_period(s)
#         return np.rad2deg(
#             np.arctan(self.tan_alpha_delta_half * np.sin(2 * np.pi * x) + self.tan_alpha_shift)
#         )
# class FourierSeriesFoldProfile(BaseFoldProfile):


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


def trigo_fold_profile(s, origin, wavelength, inflectionpointangle):
    """

    Parameters
    ----------
    s
    origin
    wavelength
    inflectionpointangle

    Returns
    -------

    """
    tan_alpha_delta_half = np.tan(inflectionpointangle)
    tan_alpha_shift = 0
    x = (s - origin) / wavelength
    return np.rad2deg(np.arctan(tan_alpha_delta_half * np.sin(2 * np.pi * x) + tan_alpha_shift))
