import logging

import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


def find_peaks_and_troughs(x, y):
    """

    Parameters
    ----------
    x np.array or list
        x axis data for plot
    y np.array or list
        y axis data for plot
    Returns
    -------
    (np.array, np.array)
    Notes
    -----
    Returns the loations of maxima/minima on the curve using finite
    difference forward/backwards
    finding the change in derivative
    """
    if len(x) != len(y):
        return False
    pairsx = []
    pairsy = []
    # #TODO numpyize
    for i in range(0, len(x)):
        if i < 1:
            pairsx.append(x[i])
            pairsy.append(y[i])

            continue
        if i > len(x) - 2:
            pairsx.append(x[i])
            pairsy.append(y[i])
            continue
        left_grad = (y[i - 1] - y[i]) / (x[i - 1] - x[i])
        right_grad = (y[i] - y[i + 1]) / (x[i] - x[i + 1])
        if np.sign(left_grad) != np.sign(right_grad):
            pairsx.append(x[i])
            pairsy.append(y[i])
    return pairsx, pairsy


class SVariogram():
    """
    The SVariogram is an experimental semi-variogram.
    """

    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        self.dist = np.abs(self.xdata[:, None] - self.xdata[None, :])
        self.variance_matrix = (self.ydata[:, None] - self.ydata[None, :]) ** 2
        self.lags = None
        self.variogram = None

    def calc_semivariogram(self, lag = None, nlag = None, lags = None):
        """
        Calculate a semi-variogram for the x and y data for this object.
        You can specify the lags as an array or specify the step size and
        number of steps.
        If neither are specified then the lags are created to be the average
        spacing of the data

        Parameters
        ----------
        step:   float
                lag distance for the s-variogram
        nstep:  int
                number of lags for the s-variogram
        lags:   array
                num

        Returns
        -------

        """
        logger.info("Calculating S-Variogram")
        if lag is not None:
            step = lag
            logger.info("Using lag: %f kwarg for S-variogram"%step)

        if nlag is not None:
            nstep = nlag
            logger.info("Using nlag %i kwarg for s-variogram"%nstep)

            self.lags = np.arange(step / 2., nstep * step, step)

        if nlag is None and lag is not None:
            nstep = int(
                np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            logger.info("Using lag kwarg but calculating nlag as %i for s-variogram"%nstep)

            self.lags = np.arange(step / 2., nstep * step, step)

        if lags is not None:
            self.lags = lags

        if self.lags is None:
            # time to guess the step size
            # find the average distance between elements in the input data
            d = np.copy(self.dist)
            d[d == 0] = np.nan

            step = np.nanmean(np.nanmin(d, axis=1))*4.
            # find number of steps to cover range in data
            nstep = int(
                np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            self.lags = np.arange(step / 2., nstep * step, step)
            logger.info("Using average minimum nearest neighbour distance "
                        "as lag distance size {} and using {} lags".format(step,nstep))
        tol = self.lags[1] - self.lags[0]
        self.variogram = np.zeros(self.lags.shape)
        self.variogram[:] = np.nan
        npairs = np.zeros(self.lags.shape)
        for i in range(len(self.lags)):
            logic = np.logical_and(self.dist > self.lags[i]
                                   - tol / 2.,
                                   self.dist < self.lags[i] + tol / 2.)
            npairs[i] = np.sum(logic.astype(int))
            if npairs[i] > 0:
                self.variogram[i] = np.mean(self.variance_matrix[logic])
        return self.lags, self.variogram, npairs

    def find_wavelengths(self, **kwargs):
        """
        Picks the wavelengths of the fold by finding the maximum and
        minimums of the s-variogram
        the fold wavelength is the first minimum but it is more reliable to
        use the first maximum
        as the estimate of the wavelength.

        Parameters
        ----------
        kwargs : object
        """
        h, var, npairs = self.calc_semivariogram(**kwargs)

        px, py = find_peaks_and_troughs(h, var)

        averagex = []
        averagey = []
        for i in range(len(px) - 1):
            averagex.append((px[i] + px[i + 1]) / 2.)
            averagey.append((py[i] + py[i + 1]) / 2.)
            i += 1  # iterate twice
        # find the extrema of the average curve
        px2, py2 = find_peaks_and_troughs(averagex, averagey)
        wl1 = 0.
        wl1py = 0.
        for i in range(len(px)):
            if i > 0 and i < len(px) - 1:
                if py[i] > 10:

                    if py[i - 1] < py[i] * .7:
                        if py[i + 1] < py[i] * .7:
                            wl1 = px[i]
                            if wl1 > 0.:
                                wl1py = py[i]
                                break
        wl2 = 0.
        for i in range(len(px2)):
            if i > 0 and i < len(px2) - 1:
                if py2[i - 1] < py2[i] * .90:
                    if py2[i + 1] < py2[i] * .90:
                        wl2 = px2[i]
                        if wl2 > 0. and wl2 > wl1 * 2 and wl1py < py2[i]:
                            break
        if wl1 == 0.0 and wl2 == 0.0:
            return 2 * (np.max(self.xdata) - np.min(self.xdata)), 0.
        if np.isclose(wl1, 0.0):
            return np.array([wl2 * 2., wl1 * 2.])
        # wavelength is 2x the peak on the curve
        return np.array([wl1 * 2., wl2 * 2.])
