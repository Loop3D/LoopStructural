from typing import List, Optional, Union

import numpy as np
import numpy.typing as npt

from .....utils import getLogger
from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile

logger = getLogger(__name__)


class TrigoFoldRotationAngleProfile(BaseFoldRotationAngleProfile):
    def __init__(
        self,
        rotation_angle: Optional[npt.NDArray[np.float64]] = None,
        fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
        origin: float = 0,
        wavelength: float = 0,
        inflectionpointangle: float = 0,
    ):
        """The fold frame function using the trigo profile from Laurent 2016

        Parameters
        ----------
        rotation_angle : npt.NDArray[np.float64], optional
            the calculated fold rotation angle from observations in degrees, by default None
        fold_frame_coordinate : npt.NDArray[np.float64], optional
            fold frame coordinate scalar field value, by default None
        origin : float, optional
            phase shift parameter, by default 0
        wavelength : float, optional
            wavelength of the profile, by default 0
        inflectionpointangle : float, optional
            height of the profile, tightness of fold, by default 0
        """
        super().__init__(rotation_angle, fold_frame_coordinate)
        self._origin = origin
        self._wavelength = wavelength
        self._inflectionpointangle = inflectionpointangle

    @property
    def origin(self):
        return self._origin

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def inflectionpointangle(self):
        return self._inflectionpointangle

    @origin.setter
    def origin(self, value):
        if np.isfinite(value):
            self.notify_observers()

            self._origin = value
        else:
            raise ValueError("origin must be a finite number")

    @wavelength.setter
    def wavelength(self, value):
        if np.isfinite(value):
            self.notify_observers()

            self._wavelength = value
        else:
            raise ValueError("wavelength must be a finite number")

    @inflectionpointangle.setter
    def inflectionpointangle(self, value):
        if np.isfinite(value):
            if value < np.deg2rad(-90) or value > np.deg2rad(90):
                logger.error(f"Inflection point angle is {np.rad2deg(value)} degrees")
                raise ValueError("inflectionpointangle must be between 0 and 90")
            self.notify_observers()
            self._inflectionpointangle = value

        else:
            raise ValueError("inflectionpointangle must be a finite number")

    @property
    def params(self):
        return {
            "origin": self.origin,
            "wavelength": self.wavelength,
            "inflectionpointangle": self.inflectionpointangle,
        }

    @staticmethod
    def _function(s, origin, wavelength, inflectionpointangle):
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
        return tan_alpha_delta_half * np.sin(2 * np.pi * x) + tan_alpha_shift

    # def __call__(self, fold_frame_coordinate):
    #     return np.rad2deg(
    #         np.tan(
    #             self._function(
    #                 fold_frame_coordinate, self.origin, self.wavelength, self.inflectionpointangle
    #             )
    #         )
    #     )

    def calculate_misfit(
        self, rotation_angle: np.ndarray, fold_frame_coordinate: np.ndarray
    ) -> np.ndarray:
        return super().calculate_misfit(rotation_angle, fold_frame_coordinate)

    def update_params(self, params: Union[List, npt.NDArray[np.float64]]) -> None:
        self.origin = params[0]
        self.wavelength = params[1]
        self.inflectionpointangle = params[2]

    def initial_guess(
        self,
        wavelength: Optional[float] = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict = {},
        reset: bool = True,
    ):
        # reset the fold paramters before fitting
        # otherwise use the current values to fit
        if reset:
            self.origin = 0
            self.wavelength = 0
            self.inflectionpointangle = np.deg2rad(45)  # otherwise there is a numerical error
        if calculate_wavelength:
            self.wavelength = self.estimate_wavelength(svariogram_parameters=svariogram_parameters)
        if wavelength is not None:
            self.wavelength = wavelength
        guess = [
            self.fold_frame_coordinate.mean(),
            self.wavelength,
            np.max(np.arctan(np.deg2rad(self.rotation_angle))),
        ]
        return np.array(guess)
