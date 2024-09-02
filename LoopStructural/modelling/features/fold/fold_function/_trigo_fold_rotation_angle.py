from _base_fold_rotation_angle import BaseFoldRotationAngle
import numpy as np

from scipy.optimize import curve_fit
from ...svariogram import SVariogram
from ....utils import getLogger


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


class TrigoFoldRotationAngle(BaseFoldRotationAngle):
    def __init__(self, rotation_angle, fold_frame_coordinate):
        self.rotation_angle = rotation_angle
        self.fold_frame_coordinate = fold_frame_coordinate
        self._origin = 0
        self._wavelength = 0
        self._inflectionpointangle = 0

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
            self._origin = value
        else:
            raise ValueError("origin must be a finite number")

    @wavelength.setter
    def wavelength(self, value):
        if np.isfinite(value):
            self._wavelength = value
        else:
            raise ValueError("wavelength must be a finite number")

    @inflectionpointangle.setter
    def inflectionpointangle(self, value):
        if np.isfinite(value):
            if value < np.deg2rad(-90) or value > np.deg2rad(90):
                raise ValueError("inflectionpointangle must be between 0 and 90")
            self._inflectionpointangle = value

        else:
            raise ValueError("inflectionpointangle must be a finite number")

    def params(self):
        return {
            "origin": self.origin,
            "wavelength": self.wavelength,
            "inflectionpointangle": self.inflectionpointangle,
        }

    def fit(self, rotation_angle, fold_frame_coordinate, params={}) -> bool:
        """Fit the trigo fold profile

        Parameters
        ----------
        rotation_angle : _type_
            _description_
        fold_frame_coordinate : _type_
            _description_
        params : dict, optional
            - wl, fold wavelength
            - origin, origin of the fold profile
            - inflectionpointangle, inflection point angle
            - svariogram_params, see :func:`~LoopStructural.modelling.fold.Svariogram.calc_semivariogram`, by default {}

        Returns
        -------
        _type_
            _description_
        """
        success = False

        wl = params.get('wl', None)
        if wl is None:
            svariogram = SVariogram(self.fold_frame_coordinate, self.rotation_angle)
            wl = svariogram.find_wavelengths(**params.get('svariogram_params', {}))
        origin = params.get('origin', 0)
        inflectionpointangle = params.get('inflectionpointangle', 45)  ## default to 45 degrees
        guess = np.array([origin, wl, np.deg2rad(inflectionpointangle)])
        logger.info(f"Guess: {guess[0]} {guess[1]} {guess[2]} ")
        # mask nans
        mask = np.logical_or(~np.isnan(self.fold_frame_coordinate), ~np.isnan(self.rotation_angle))
        logger.info(
            f"There are {np.sum(~mask)} nans for the fold limb rotation angle and { np.sum(mask)} observations"
        )
        if np.sum(mask) < len(guess):
            logger.error(
                "Not enough data points to fit Fourier series setting " "fold rotation angle" "to 0"
            )
            self.fold_rotation_function = lambda x: np.zeros(x.shape)
        else:
            try:
                # try fitting using wavelength guess
                popt, pcov = curve_fit(
                    trigo_fold_profile,
                    self.fold_frame_coordinate[mask],
                    np.tan(np.deg2rad(self.rotation_angle[mask])),
                    guess,
                )
                success = True
            except RuntimeError:
                try:
                    # if fitting failed, try with just 0s
                    logger.info("Running curve fit without initial guess")
                    popt, pcov = curve_fit(
                        trigo_fold_profile,
                        self.fold_frame_coordinate[mask],
                        np.tan(np.deg2rad(self.rotation_angle[mask])),
                    )
                    success = True
                except RuntimeError:
                    # otherwise set the fourier series parameters to 0
                    popt = guess
                    success = False
                    logger.error("Could not fit curve to S-Plot, check the wavelength")
            self._origin = popt[0]
            self._wavelength = popt[1]
            self._inflectionpointangle = popt[2]
            return success

    def __call__(self, fold_frame_coordinate):
        return trigo_fold_profile(
            fold_frame_coordinate, self.origin, self.wavelength, self.inflectionpointangle
        )
