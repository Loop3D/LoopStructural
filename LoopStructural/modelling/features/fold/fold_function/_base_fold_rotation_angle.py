from abc import ABCMeta, abstractmethod
from ast import List
from typing import Union, Optional
import numpy as np
import numpy.typing as npt
from .._svariogram import SVariogram
from scipy.optimize import curve_fit

from .....utils import getLogger

logger = getLogger(__name__)


class BaseFoldRotationAngleProfile(metaclass=ABCMeta):

    def __init__(
        self,
        rotation_angle: Optional[npt.NDArray[np.float64]] = None,
        fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
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

    @abstractmethod
    def fit(self, params: dict = {}) -> bool:
        """Method is called to setup the fold rotation angle for the given data
        any parameters that are required to fit the fold rotation angle can be injected
        using the params dictionary

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
        success : bool
            returns True if the fit was successful
        """
        pass

    # @abstractmethod
    # def __call__(self, fold_frame_coordinate: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    #     """This is called during the fold interpolation to return the fold rotation angle for a given
    #     fold frame coordinate.
    #     Should be able to be called with a single value or an array of values

    #     Parameters
    #     ----------
    #     fold_frame_coordinate : Union[float, np.ndarray]
    #         fold frame coordinate to return a rotation angle

    #     Returns
    #     -------
    #     Union[float, np.ndarray]
    #         fold rotation angle
    #     """
    #     pass

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
        self, svariogram_parameters: dict = {}, wavelength_number: int = 1
    ) -> Union[float, np.ndarray]:
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
        svariogram = SVariogram(self.fold_frame_coordinate, self.rotation_angle)
        wl = svariogram.find_wavelengths(**svariogram_parameters)
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

    def fit(self, params: dict = {}) -> bool:
        """

        Parameters
        ----------
        params : dict, optional
            _description_, by default {}

        Returns
        -------
        bool
            _description_
        """
        success = False
        if self.rotation_angle is None or self.fold_frame_coordinate is None:
            logger.error("Fold rotation angle and fold frame coordinate must be set")
            return False
        guess = params.get(
            "guess",
            self.initial_guess(
                params.get("reset", False),
                params.get("svariogram_parameters", {}),
                params.get("calculate_wavelength", True),
            ),
        )
        mask = np.logical_or(~np.isnan(self.fold_frame_coordinate), ~np.isnan(self.rotation_angle))
        logger.info(f"Percentage of points not used {np.sum(~mask)/len(mask)*100}")
        try:
            logger.info(f"Trying to fit fold rotation angle with guess {guess}")
            logger.info(f"Fold profile type: {self.__class__.__name__}")
            res = curve_fit(
                self._function,
                self.fold_frame_coordinate[mask],
                np.arctan(np.deg2rad(self.rotation_angle[mask])),
                p0=guess,
                full_output=True,
            )
            logger.info(f'Fit results: {res[0]}')
            guess = res[0]
            logger.info(res[3])
            success = True
        except Exception as e:
            logger.error("Could not fit curve to S-Plot, check the wavelength")

        self.update_params(guess)
        return success

    @abstractmethod
    def update_params(self, params: Union[List, npt.NDArray[np.float64]]) -> None:
        """Update the parameters of the fold rotation angle function

        Parameters
        ----------
        params : dict
            parameters to update
        """
        pass

    @abstractmethod
    def initial_guess(
        self,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict = {},
        reset: bool = False,
    ) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        selfcalculate_wavelength : bool, optional
            _description_, by default True
        svariogram_parameters : dict, optional
            _description_, by default {}
        reset : bool, optional
            _description_, by default False

        Returns
        -------
        np.ndarray
            _description_
        """
        pass

    @staticmethod
    @abstractmethod
    def _function(s, *args, **kwargs):
        """This is the function that is used to calculate the fold rotation angle
        for a given fold frame coordinate
        it is not called directly but is used by the __call__ method

        Parameters
        ----------
        s
        *args

        Returns
        -------
        _description_
        """
        pass

    def plot(self, ax=None, **kwargs):
        """Plot the fold rotation angle function

        Parameters
        ----------
        ax : _description_, optional
            _description_, by default None
        **kwargs
            passed to matplotlib plot
        """
        if ax is None:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
        ax.plot(self.evaluation_points, self(self.evaluation_points), **kwargs)
        return ax

    def __call__(self, s):
        return np.rad2deg(np.tan(self._function(s, **self.params)))
