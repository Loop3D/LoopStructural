"""Semi-variogram for estimating fold wavelengths from orientation data."""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from loop_common.logging import get_logger

logger = get_logger(__name__)


def find_peaks_and_troughs(x: np.ndarray, y: np.ndarray) -> Tuple[List, List]:
    """Return x/y positions of local maxima and minima using finite differences."""
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    pairsx: List = []
    pairsy: List = []
    for i in range(len(x)):
        if i < 1 or i > len(x) - 2:
            if not np.isnan(y[i]):
                pairsx.append(x[i])
                pairsy.append(y[i])
            continue
        if np.isnan(y[i - 1]) or np.isnan(y[i]) or np.isnan(y[i + 1]):
            continue
        left_grad = (y[i - 1] - y[i]) / (x[i - 1] - x[i])
        right_grad = (y[i] - y[i + 1]) / (x[i] - x[i + 1])
        if np.sign(left_grad) != np.sign(right_grad):
            pairsx.append(x[i])
            pairsy.append(y[i])
    return pairsx, pairsy


class SVariogram:
    """Experimental semi-variogram used to estimate fold wavelengths."""

    def __init__(self, xdata: np.ndarray, ydata: np.ndarray):
        self.xdata = np.asarray(xdata)
        self.ydata = np.asarray(ydata)
        mask = np.logical_or(np.isnan(self.xdata), np.isnan(self.ydata))
        self.xdata = self.xdata[~mask]
        self.ydata = self.ydata[~mask]
        self.dist = np.abs(self.xdata[:, None] - self.xdata[None, :])
        self.variance_matrix = (self.ydata[:, None] - self.ydata[None, :]) ** 2
        self.lags: Optional[np.ndarray] = None
        self.variogram: Optional[np.ndarray] = None
        self.wavelength_guesses: List = []

    def initialise_lags(self, step: Optional[float] = None, nsteps: Optional[int] = None):
        if nsteps is not None and step is not None:
            self.lags = np.arange(step / 2.0, nsteps * step, step)
        elif step is not None:
            nsteps = int(np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            self.lags = np.arange(step / 2.0, nsteps * step, step)

        if self.lags is None:
            d = np.copy(self.dist)
            d[d == 0] = np.nan
            step = np.nanmean(np.nanmin(d, axis=1)) * 4.0
            nsteps = int(np.ceil((np.nanmax(self.xdata) - np.nanmin(self.xdata)) / step))
            if nsteps > 200:
                logger.warning(f"Variogram has too many steps: {nsteps}, capping at 200")
                maximum = step * nsteps
                nsteps = 200
                step = maximum / nsteps
            self.lags = np.arange(step / 2.0, nsteps * step, step)

    def calc_semivariogram(
        self,
        step: Optional[float] = None,
        nsteps: Optional[int] = None,
        lags: Optional[np.ndarray] = None,
    ):
        if lags is not None:
            self.lags = lags
        self.initialise_lags(step, nsteps)
        if self.lags is None:
            raise ValueError(
                "Cannot determine variogram step size; specify step or nsteps."
            )
        tol = self.lags[1] - self.lags[0]
        self.variogram = np.full(self.lags.shape, np.nan)
        npairs = np.zeros(self.lags.shape)
        for i in range(len(self.lags)):
            logic = (self.dist > self.lags[i] - tol / 2.0) & (
                self.dist < self.lags[i] + tol / 2.0
            )
            npairs[i] = np.sum(logic)
            if npairs[i] > 0:
                self.variogram[i] = np.mean(self.variance_matrix[logic])
        return self.lags, self.variogram, npairs

    def find_wavelengths(
        self,
        step: Optional[float] = None,
        nsteps: Optional[int] = None,
        lags: Optional[np.ndarray] = None,
    ) -> List:
        h, var, _npairs = self.calc_semivariogram(step=step, nsteps=nsteps, lags=lags)
        px, py = find_peaks_and_troughs(h, var)

        averagex: List = []
        averagey: List = []
        for i in range(len(px) - 1):
            averagex.append((px[i] + px[i + 1]) / 2.0)
            averagey.append((py[i] + py[i + 1]) / 2.0)

        res = find_peaks_and_troughs(np.array(averagex), np.array(averagey))
        px2, py2 = res

        wl1 = 0.0
        wl1py = 0.0
        for i in range(len(px)):
            if 0 < i < len(px) - 1 and py[i] > 10:
                if py[i - 1] < py[i] * 0.7 and py[i + 1] < py[i] * 0.7:
                    wl1 = px[i]
                    if wl1 > 0.0:
                        wl1py = py[i]
                        break

        wl2 = 0.0
        for i in range(len(px2)):
            if 0 < i < len(px2) - 1:
                if py2[i - 1] < py2[i] * 0.90 and py2[i + 1] < py2[i] * 0.90:
                    wl2 = px2[i]
                    if wl2 > 0.0 and wl2 > wl1 * 2 and wl1py < py2[i]:
                        break

        if wl1 == 0.0 and wl2 == 0.0:
            logger.warning("Could not auto-estimate wavelength; using 2× data range")
            self.wavelength_guesses = [2 * (np.max(self.xdata) - np.min(self.xdata)), 0.0]
            return self.wavelength_guesses
        if np.isclose(wl1, 0.0):
            self.wavelength_guesses = [wl2 * 2.0, 0.0]
            return [wl2 * 2.0]
        self.wavelength_guesses = [wl1 * 2.0, wl2 * 2.0]
        return self.wavelength_guesses
