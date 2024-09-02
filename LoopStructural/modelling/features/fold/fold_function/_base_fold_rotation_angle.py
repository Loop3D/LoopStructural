from abc import ABCMeta, abstractmethod
from typing import Union
import numpy as np


class BaseFoldRotationAngle(ABCMeta):
    def __init__(self):
        """Base class for fold rotation angle functions."""

    @abstractmethod
    def fit(
        self, rotation_angle: np.ndarray, fold_frame_coordinate: np.ndarray, params: dict = {}
    ) -> bool:
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

    @abstractmethod
    def __call__(self, fold_frame_coordinate: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """This is called during the fold interpolation to return the fold rotation angle for a given
        fold frame coordinate.
        Should be able to be called with a single value or an array of values

        Parameters
        ----------
        fold_frame_coordinate : Union[float, np.ndarray]
            fold frame coordinate to return a rotation angle

        Returns
        -------
        Union[float, np.ndarray]
            fold rotation angle
        """
        pass

    @abstractmethod
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

        pass

    
