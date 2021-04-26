import logging

import numpy as np
from scipy.optimize import curve_fit

from LoopStructural.modelling.fold.fold_rotation_angle_feature import \
    fourier_series
from LoopStructural.modelling.fold import SVariogram

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class FoldRotationAngle:
    """

    """
    def __init__(self, rotation_angle, fold_frame_coordinate, svario=False):
        """

        Parameters
        ----------
        rotation_angle
        fold_frame_coordinate
        svario
        """
        self.rotation_angle = rotation_angle
        self.fold_frame_coordinate = fold_frame_coordinate
        self.fold_rotation_function = None
        self.svario = None
        if svario:
            self.svario = SVariogram(self.fold_frame_coordinate, self.rotation_angle)
        self.fitted_params = None

    def fit_fourier_series(self, wl = None, lags = None, nlag = None, lag = None, skip_variogram=False,**kwargs):
        """

        Parameters
        ----------
        wl
        lags
        nlag
        lag

        Returns
        -------

        """
        if self.svario is None:
            self.svario = SVariogram(self.fold_frame_coordinate,
                                     self.rotation_angle)
        if skip_variogram == False:
            self.svario.calc_semivariogram(lags=lags, nlag=nlag, lag=lag)
        if wl is None:
            wl = self.svario.find_wavelengths(lags=lags, nlag=nlag, lag=lag)
            # for now only consider single fold wavelength
            wl = wl[0]
        guess = np.zeros(4)
        guess[3] = wl  # np.max(limb_wl)
        logger.info(
            'Guess: %f %f %f %f' % (guess[0], guess[1], guess[2], guess[3]))
        # mask nans
        mask = np.logical_or(~np.isnan(self.fold_frame_coordinate), ~np.isnan(self.rotation_angle))
        logger.info(
            "There are %i nans for the fold limb rotation angle and "
            "%i observations" % (np.sum(~mask), np.sum(mask)))
        if np.sum(mask) < len(guess):
            logger.error(
                "Not enough data points to fit Fourier series setting "
                "fold rotation angle"
                "to 0")
            self.fold_rotation_function = lambda x: np.zeros(x.shape)
        else:
            try:
                # try fitting using wavelength guess
                popt, pcov = curve_fit(fourier_series,
                                       self.fold_frame_coordinate[mask],
                                       np.tan(np.deg2rad(self.rotation_angle[mask])),
                                       guess)
            except RuntimeError:
                try:
                    # if fitting failed, try with just 0s
                    logger.info("Running curve fit without initial guess")
                    popt, pcov = curve_fit(fourier_series,
                                           self.fold_frame_coordinate[mask],
                                           np.tan(np.deg2rad(self.rotation_angle[mask])))
                except RuntimeError:
                    # otherwise set the fourier series parameters to 0
                    popt = guess
                    logger.error(
                        "Could not fit curve to S-Plot, check the wavelength")
            logger.info(
                'Fitted: %f %f %f %f' % (popt[0], popt[1], popt[2], popt[3]))
            self.fold_rotation_function = lambda x: np.rad2deg(
                np.arctan(
                    fourier_series(x, popt[0], popt[1], popt[2], popt[3])))
            self.fitted_params = popt

    def __call__(self, fold_frame_coordinate):
        """

        Parameters
        ----------
        fold_frame_coordinate

        Returns
        -------

        """
        return self.fold_rotation_function(fold_frame_coordinate)

    def calculate_misfit(self):
        """

        Returns
        -------

        """
        return np.tan(np.deg2rad(self.rotation_angle)) - np.tan(np.deg2rad(
            self.__call__(self.fold_frame_coordinate)))

    def set_function(self, function):
        """

        Parameters
        ----------
        function

        Returns
        -------

        """
        self. fold_rotation_function = function

    def find_hinges(self,range,step):
        
        import scipy.optimize as optimize
        def fra(x):
            x = np.array([x])
            return self.__call__(x)
        roots = []
        x = range[0]
        while x < range[1]:
            result = optimize.root_scalar(fra,bracket=[x,x+step])
            roots.append(result.root)
            x+=step
        return roots