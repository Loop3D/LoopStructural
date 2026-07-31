from __future__ import annotations

from abc import ABCMeta, abstractmethod
from ast import List

import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit

from .....utils import getLogger
from .._svariogram import SVariogram

logger = getLogger(__name__)


class BaseFoldRotationAngleProfile(metaclass=ABCMeta):
    def __init__(
        self,
        rotation_angle: npt.NDArray[np.float64] | None = None,
        fold_frame_coordinate: npt.NDArray[np.float64] | None = None,
    ):
        """Base class for fold rotation angle functions

        Parameters
        ----------
        rotation_angle : npt.NDArray[np.float64], optional
            the calculated fold rotation angle from observations in degrees, by default None
        fold_frame_coordinate : npt.NDArray[np.float64], optional
            fold frame coordinate scalar field value, by default None
        """

        self.rotation_angle = rotation_angle
        self.fold_frame_coordinate = fold_frame_coordinate
        self._evaluation_points = None
        self._observers = []
        self._svariogram = None

    @property
    def svario(self) -> SVariogram:
        if self.fold_frame_coordinate is None or self.rotation_angle is None:
            raise ValueError("Fold rotation angle and fold frame coordinate must be set")
        if self._svariogram is None:
            self._svariogram = SVariogram(self.fold_frame_coordinate, self.rotation_angle)
        return self._svariogram

    @svario.setter
    def svario(self, value: SVariogram):
        if isinstance(value, SVariogram):
            self._svariogram = value
        else:
            logger.error("svario must be an instance of SVariogram")
            raise ValueError("svario must be an instance of SVariogram")

    def add_observer(self, watcher):
        self._observers.append(watcher)

    def notify_observers(self):
        for observer in self._observers:
            observer.set_not_up_to_date(self)

    def calculate_misfit(
        self,
        rotation_angle: np.ndarray,
        fold_frame_coordinate: np.ndarray,
    ) -> np.ndarray:
        """Calculate the rotation angle for the fold frame coordinate and return the misfit

        Parameters
        ----------
        params : dict, optional
            Any parameters required to fit the function to the data, by default {}
        rotation_angle : np.ndarray
            fold rotation angle in degrees
        fold_frame_coordinate : np.ndarray
            fold frame coordinate

        Returns
        -------
        misfit : np.ndarray
            returns misfit in degrees"""
        return np.tan(np.deg2rad(rotation_angle)) - np.tan(
            np.deg2rad(self.__call__(fold_frame_coordinate))
        )

    def estimate_wavelength(
        self, svariogram_parameters: dict | None = None, wavelength_number: int = 1
    ) -> float | np.ndarray:
        """Estimate the wavelength of the fold profile using the svariogram parameters

        Parameters
        ----------
        svariogram_parameters : dict
            svariogram parameters

        Returns
        -------
        float
            estimated wavelength
        """
        if svariogram_parameters is None:
            svariogram_parameters = {}
        wl = self.svario.find_wavelengths(**svariogram_parameters)
        if wavelength_number == 1:
            return wl[0]
        return wl

    @property
    def evaluation_points(self):
        """Return the evaluation points for the fold rotation angle function

        Returns
        -------
        np.ndarray
            evaluation points
        """
        if self._evaluation_points is not None:
            return self._evaluation_points
        return np.linspace(
            np.min(self.fold_frame_coordinate), np.max(self.fold_frame_coordinate), 300
        )

    @evaluation_points.setter
    def evaluation_points(self, value):
        self._evaluation_points = value

    def fit(self, params: dict | None = None) -> bool:
        """Fit the fold rotation angle function to the rotation angle and fold frame
        coordinate observations using scipy curve_fit

        Parameters
        ----------
        params : dict, optional
            fitting parameters, may contain "guess", "wavelength", "reset",
            "svariogram_parameters" and "calculate_wavelength" keys used to build the
            initial guess, by default {}

        Returns
        -------
        bool
            True if the curve was successfully fit, False otherwise
        """
        if params is None:
            params = {}
        if len(self.params) > 0:
            success = False
            if self.rotation_angle is None or self.fold_frame_coordinate is None:
                logger.error("Fold rotation angle and fold frame coordinate must be set")
                return False
            guess = params.get(
                "guess",
                self.initial_guess(
                    wavelength=params.get("wavelength", None),
                    reset=params.get("reset", False),
                    svariogram_parameters=params.get("svariogram_parameters", {}),
                    calculate_wavelength=params.get("calculate_wavelength", True),
                ),
            )
            mask = np.logical_or(
                ~np.isnan(self.fold_frame_coordinate), ~np.isnan(self.rotation_angle)
            )
            logger.info(f"Percentage of points not used {np.sum(~mask)/len(mask)*100}")
            try:
                logger.info(f"Trying to fit fold rotation angle with guess {guess}")
                logger.info(f"Fold profile type: {self.__class__.__name__}")
                res = curve_fit(
                    self._function,
                    self.fold_frame_coordinate[mask],
                    np.tan(np.deg2rad(self.rotation_angle[mask])),
                    p0=guess,
                    full_output=True,
                )
                logger.info(f'Fit results: {res[0]}')
                guess = res[0]
                logger.info(res[3])
                success = True
            except (RuntimeError, ValueError, TypeError) as _e:
                logger.error(f"Could not fit curve to S-Plot, check the wavelength: {_e}")
            try:
                self.update_params(guess)
            except (ValueError, TypeError, IndexError) as _e:
                logger.error(f"Could not update parameters: {_e}")
                return False
            return success
        return True

    @abstractmethod
    def update_params(self, params: List | npt.NDArray[np.float64]) -> None:
        """Update the parameters of the fold rotation angle function

        Parameters
        ----------
        params : dict
            parameters to update
        """

    @abstractmethod
    def initial_guess(
        self,
        wavelength: float | None = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict | None = None,
        reset: bool = False,
    ) -> np.ndarray:
        """Calculate an initial guess for the parameters of the fold rotation angle function,
        optionally using the wavelength estimated from the svariogram

        Parameters
        ----------
        wavelength : float, optional
            wavelength to use for the initial guess, if None it is estimated from the
            svariogram when calculate_wavelength is True, by default None
        calculate_wavelength : bool, optional
            whether to estimate the wavelength from the svariogram, by default True
        svariogram_parameters : dict, optional
            parameters passed to the svariogram when estimating the wavelength, by default {}
        reset : bool, optional
            whether to reset any previously fitted parameters before guessing, by default False

        Returns
        -------
        np.ndarray
            initial guess of the parameters for the fold rotation angle function
        """
        if svariogram_parameters is None:
            svariogram_parameters = {}

    @staticmethod
    @abstractmethod
    def _function(s, *args, **kwargs):
        """This is the function that is used to calculate the fold rotation angle
        for a given fold frame coordinate
        it is not called directly but is used by the __call__ method

        Parameters
        ----------
        s
            fold frame coordinate to evaluate the function at
        *args
            parameters of the fold rotation angle function

        Returns
        -------
        np.ndarray
            tan of the fold rotation angle in radians at s
        """

    def plot(self, ax=None, show_data=True, **kwargs):
        """Plot the fold rotation angle function

        Parameters
        ----------
        ax : matplotlib axes, optional
            the axes to plot onto, a new figure and axes are created if None, by default None
        **kwargs
            passed to matplotlib plot
        """
        if ax is None:
            import matplotlib.pyplot as plt

            _fig, ax = plt.subplots()
        if show_data:
            ax.scatter(self.fold_frame_coordinate, self.rotation_angle, c="r")
        ax.plot(self.evaluation_points, self(self.evaluation_points), **kwargs)
        return ax

    def __call__(self, s):
        return np.rad2deg(np.arctan(self._function(s, **self.params)))
