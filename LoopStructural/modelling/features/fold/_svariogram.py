from __future__ import annotations

import numpy as np

from ....utils import getLogger

logger = getLogger(__name__)


def find_peaks_and_troughs(x: np.ndarray, y: np.ndarray) -> tuple[list, list]:
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
        raise ValueError("Cannot guess wavelength, x and y must be the same length")
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    if n == 0:
        return [], []
    # always keep the first and last point; keep interior points where the
    # sign of the finite-difference gradient changes (local max/min)
    mask = np.zeros(n, dtype=bool)
    mask[0] = True
    mask[-1] = True
    if n > 2:
        left_grad = (y[:-2] - y[1:-1]) / (x[:-2] - x[1:-1])
        right_grad = (y[1:-1] - y[2:]) / (x[1:-1] - x[2:])
        mask[1:-1] = np.sign(left_grad) != np.sign(right_grad)
    return list(x[mask]), list(y[mask])


class SVariogram:
    """
    The SVariogram is an experimental semi-variogram.
    """

    def __init__(self, xdata: np.ndarray, ydata: np.ndarray):
        self.xdata = np.asarray(xdata)
        self.ydata = np.asarray(ydata)
        mask = np.logical_or(np.isnan(self.xdata), np.isnan(self.ydata))
        self.xdata = self.xdata[~mask]
        self.ydata = self.ydata[~mask]
        ## maybe check that xdata is not too big here, if it is then this should be done
        ## using another library maybe gstools?
        self.dist = np.abs(self.xdata[:, None] - self.xdata[None, :])
        self.variance_matrix = (self.ydata[:, None] - self.ydata[None, :]) ** 2
        self.lags = None
        self.variogram = None
        self.wavelength_guesses = []

    def initialise_lags(self, step: float | None = None, nsteps: int | None = None):
        """
        Initialise the lags for the s-variogram

        Parameters
        ----------
        lag: float
            lag distance for the s-variogram
        nlag: int
            number of lags for the s-variogram

        Returns
        -------

        """

        if nsteps is not None and step is not None:
            logger.info(f"Using nlag {nsteps} kwarg for s-variogram")
            logger.info(f"Using lag: {step} kwarg for S-variogram")
            self.lags = np.arange(step / 2.0, nsteps * step, step)

        if nsteps is None and step is not None:
            nsteps = int(np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            logger.info(f"Using lag kwarg but calculating nlag as {nsteps} for s-variogram")

            self.lags = np.arange(step / 2.0, nsteps * step, step)

        if self.lags is None:
            # time to guess the step size
            # find the average distance between elements in the input data
            d = np.copy(self.dist)
            d[d == 0] = np.nan

            step = np.nanmean(np.nanmin(d, axis=1)) * 4.0
            # find number of steps to cover range in data
            nsteps = int(np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            if nsteps > 200:
                logger.warning(f"Variogram has too many steps: {nsteps}, using 200")
                maximum = step * nsteps
                nstep = 200
                step = maximum / nstep
            self.lags = np.arange(step / 2.0, nsteps * step, step)
            logger.info(
                f"Using average minimum nearest neighbour distance as lag distance size {step} and using {nsteps} lags"
            )

    def calc_semivariogram(
        self,
        step: float | None = None,
        nsteps: int | None = None,
        lags: np.ndarray | None = None,
    ):
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
        if lags is not None:
            self.lags = lags
        self.initialise_lags(step, nsteps)
        if self.lags is None:
            raise ValueError(
                "S-Variogram cannot calculate the variogram step size, please specify either step or nsteps"
            )
        tol = self.lags[1] - self.lags[0]
        self.variogram = np.zeros(self.lags.shape)
        self.variogram[:] = np.nan
        npairs = np.zeros(self.lags.shape)
        for i in range(len(self.lags)):
            logic = np.logical_and(
                self.dist > self.lags[i] - tol / 2.0,
                self.dist < self.lags[i] + tol / 2.0,
            )
            npairs[i] = np.sum(logic.astype(int))
            if npairs[i] > 0:
                self.variogram[i] = np.mean(self.variance_matrix[logic])
        return self.lags, self.variogram, npairs

    def find_wavelengths(
        self,
        step: float | None = None,
        nsteps: int | None = None,
        lags: np.ndarray | None = None,
    ) -> list:
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
        h, var, _npairs = self.calc_semivariogram(step=step, nsteps=nsteps, lags=lags)

        px, py = find_peaks_and_troughs(h, var)

        px_arr = np.asarray(px)
        py_arr = np.asarray(py)
        if len(px_arr) > 1:
            averagex = list((px_arr[:-1] + px_arr[1:]) / 2.0)
            averagey = list((py_arr[:-1] + py_arr[1:]) / 2.0)
        else:
            averagex = []
            averagey = []
        # find the extrema of the average curve
        res = find_peaks_and_troughs(np.array(averagex), np.array(averagey))
        px2, py2 = res
        logger.info(f"Found {len(px2)} peaks and troughs in the s-variogram")
        for i in range(len(px2)):
            logger.info(f"Peak {i}: {px2[i]} {py2[i]}")
        wl1 = 0.0
        wl1py = 0.0
        for i in range(len(px)):
            if i > 0 and i < len(px) - 1 and py[i] > 10 and py[i - 1] < py[i] * 0.7:
                if py[i + 1] < py[i] * 0.7:
                    wl1 = px[i]
                    if wl1 > 0.0:
                        wl1py = py[i]
                        break
        wl2 = 0.0
        for i in range(len(px2)):
            if i > 0 and i < len(px2) - 1 and py2[i - 1] < py2[i] * 0.90:
                if py2[i + 1] < py2[i] * 0.90:
                    wl2 = px2[i]
                    if wl2 > 0.0 and wl2 > wl1 * 2 and wl1py < py2[i]:
                        break
        if wl1 == 0.0 and wl2 == 0.0:
            logger.warning(
                'Could not automatically guess the wavelength, using 2x the range of the data'
            )
            self.wavelength_guess = [2 * (np.max(self.xdata) - np.min(self.xdata)), 0.0]
            return self.wavelength_guess
        if np.isclose(wl1, 0.0):
            self.wavelength_guess = np.array([wl2 * 2.0, wl1 * 2.0])
            return [wl2]
        # wavelength is 2x the peak on the curve
        self.wavelength_guess = [wl1 * 2.0, wl2 * 2.0]
        return self.wavelength_guess
