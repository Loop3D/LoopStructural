from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .....utils import getLogger
from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile

logger = getLogger(__name__)


class FourierSeriesFoldRotationAngleProfile(BaseFoldRotationAngleProfile):
    def __init__(
        self,
        rotation_angle: npt.NDArray[np.float64] | None = None,
        fold_frame_coordinate: npt.NDArray[np.float64] | None = None,
        c0=0,
        c1=0,
        c2=0,
        w=1,
    ):
        """Fold rotation angle profile defined by a truncated Fourier series
        c0 + c1*cos(2*pi/w * x) + c2*sin(2*pi/w * x)

        Parameters
        ----------
        rotation_angle : Optional[npt.NDArray[np.float64]], optional
            the calculated fold rotation angle from observations in degrees, by default None
        fold_frame_coordinate : Optional[npt.NDArray[np.float64]], optional
            fold frame coordinate scalar field value, by default None
        c0 : int, optional
            mean value coefficient of the Fourier series, by default 0
        c1 : int, optional
            cosine coefficient of the Fourier series, by default 0
        c2 : int, optional
            sine coefficient of the Fourier series, by default 0
        w : int, optional
            wavelength of the Fourier series, by default 1
        """
        super().__init__(rotation_angle, fold_frame_coordinate)
        self._c0 = c0
        self._c1 = c1
        self._c2 = c2
        self._w = w

    @property
    def c0(self):
        return self._c0

    @c0.setter
    def c0(self, value):
        self.notify_observers()
        self._c0 = value

    @property
    def c1(self):
        return self._c1

    @c1.setter
    def c1(self, value):
        self.notify_observers()
        self._c1 = value

    @property
    def c2(self):
        return self._c2

    @c2.setter
    def c2(self, value):
        self.notify_observers()
        self._c2 = value

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        if value <= 0:
            raise ValueError('wavelength must be greater than 0')
        self.notify_observers()
        self._w = value

    @staticmethod
    def _function(x, c0, c1, c2, w):
        """Evaluate the Fourier series fold rotation angle function

        Parameters
        ----------
        x
            fold frame coordinate to evaluate the function at
        c0
            mean value coefficient of the Fourier series
        c1
            cosine coefficient of the Fourier series
        c2
            sine coefficient of the Fourier series
        w
            wavelength of the Fourier series

        Returns
        -------
        np.ndarray
            value of the Fourier series at x
        """
        v = np.array(x.astype(float))
        # v.fill(c0)
        v = c0 + c1 * np.cos(2 * np.pi / w * x) + c2 * np.sin(2 * np.pi / w * x)
        return v

    def initial_guess(
        self,
        wavelength: float | None = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict | None = None,
        reset: bool = False,
    ):
        if svariogram_parameters is None:
            svariogram_parameters = {}
        # reset the fold paramters before fitting
        # otherwise use the current values to fit
        if reset:
            self.c0 = np.mean(np.arctan(np.deg2rad(self.rotation_angle)))
            self.c1 = 0
            self.c2 = np.max(np.arctan(np.deg2rad(self.rotation_angle)))
            self.w = 1
        if calculate_wavelength:
            self.w = self.estimate_wavelength(svariogram_parameters=svariogram_parameters)
        if wavelength is not None:
            self.w = wavelength
        guess = [self.c0, self.c1, self.c2, self.w]
        return guess

    @property
    def params(self):
        return {
            "c0": self.c0,
            "c1": self.c1,
            "c2": self.c2,
            "w": self.w,
        }

    @params.setter
    def params(self, params):
        for key in params:
            if key == 'w' and params[key] <= 0:
                raise ValueError('wavelength must be greater than 0')
            setattr(self, key, params[key])
        self.c0 = params["c0"]
        self.c1 = params["c1"]
        self.c2 = params["c2"]
        self.w = params["w"]

    def update_params(self, params: list[float] | npt.NDArray[np.float64]):
        if len(params) != 4:
            raise ValueError('params must have 4 elements')
        self.c0 = params[0]
        self.c1 = params[1]
        self.c2 = params[2]
        self.w = params[3]

    def calculate_misfit(self, s, rotation_angle):
        return np.tan(np.deg2rad(rotation_angle)) - np.tan(np.deg2rad(self.__call__(s)))
