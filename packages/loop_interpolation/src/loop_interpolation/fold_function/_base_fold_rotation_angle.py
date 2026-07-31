"""Abstract base class for fold rotation-angle profiles."""
from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import List, Optional, Union

import numpy as np
import numpy.typing as npt
from loop_common.logging import get_logger
from scipy.optimize import curve_fit

from .._svariogram import SVariogram

logger = get_logger(__name__)


class BaseFoldRotationAngleProfile(metaclass=ABCMeta):
    def __init__(
        self,
        rotation_angle: Optional[npt.NDArray[np.float64]] = None,
        fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
    ):
        """Base class for callable fold-rotation-angle functions.

        Parameters
        ----------
        rotation_angle
            Observed fold rotation angles in degrees.
        fold_frame_coordinate
            Fold-frame scalar-field values at the observation locations.
        """
        self.rotation_angle = rotation_angle
        self.fold_frame_coordinate = fold_frame_coordinate
        self._evaluation_points: Optional[np.ndarray] = None
        self._observers: List = []
        self._svariogram: Optional[SVariogram] = None

    @property
    def svario(self) -> SVariogram:
        if self.fold_frame_coordinate is None or self.rotation_angle is None:
            raise ValueError("rotation_angle and fold_frame_coordinate must be set first")
        if self._svariogram is None:
            self._svariogram = SVariogram(self.fold_frame_coordinate, self.rotation_angle)
        return self._svariogram

    @svario.setter
    def svario(self, value: SVariogram):
        if not isinstance(value, SVariogram):
            raise ValueError("svario must be a SVariogram instance")
        self._svariogram = value

    def add_observer(self, watcher) -> None:
        self._observers.append(watcher)

    def notify_observers(self) -> None:
        for observer in self._observers:
            observer.set_not_up_to_date(self)

    @property
    def evaluation_points(self) -> np.ndarray:
        if self._evaluation_points is not None:
            return self._evaluation_points
        return np.linspace(
            np.min(self.fold_frame_coordinate), np.max(self.fold_frame_coordinate), 300
        )

    @evaluation_points.setter
    def evaluation_points(self, value: np.ndarray) -> None:
        self._evaluation_points = value

    def estimate_wavelength(
        self, svariogram_parameters: dict = {}, wavelength_number: int = 1
    ) -> Union[float, np.ndarray]:
        wl = self.svario.find_wavelengths(**svariogram_parameters)
        logger.info(f"Estimated fold rotation wavelength(s): {wl}")
        return wl[0] if wavelength_number == 1 else wl

    def calculate_misfit(
        self,
        rotation_angle: np.ndarray,
        fold_frame_coordinate: np.ndarray,
    ) -> np.ndarray:
        return np.tan(np.deg2rad(rotation_angle)) - np.tan(
            np.deg2rad(self.__call__(fold_frame_coordinate))
        )

    def fit(self, params: dict = {}) -> bool:
        if len(self.params) > 0:
            if self.rotation_angle is None or self.fold_frame_coordinate is None:
                logger.error("rotation_angle and fold_frame_coordinate must be set before fitting")
                return False

            # Exclude NaNs and extreme angles (|alpha| >= 89°) whose tan blows up
            # and causes the optimizer to diverge.
            mask = (
                ~np.isnan(self.fold_frame_coordinate)
                & ~np.isnan(self.rotation_angle)
                & (np.abs(self.rotation_angle) < 89.0)
            )
            n_nan = np.sum(np.isnan(self.fold_frame_coordinate) | np.isnan(self.rotation_angle))
            n_extreme = np.sum(np.abs(self.rotation_angle) >= 89.0) - np.sum(np.isnan(self.rotation_angle))
            logger.info(
                f"Fitting fold rotation angle; excluded {n_nan} NaN, {n_extreme} extreme (>=89 deg)"
            )
            x = self.fold_frame_coordinate[mask]

            # Divide the fold-frame coordinate by its range so curve_fit and the
            # variogram both operate at O(1) scale.  Only division is used (no
            # shift), so the wavelength transforms as w_orig = w_scaled * x_scale
            # while c0/c1/c2 are identical in both spaces.
            x_scale = float(np.max(x) - np.min(x)) if len(x) > 1 else 1.0
            if x_scale < 1e-12:
                x_scale = 1.0
            x_scaled = x / x_scale

            # Temporarily substitute scaled coordinates so the variogram and
            # initial_guess operate in the normalised range.
            _orig_ffc = self.fold_frame_coordinate
            _orig_svario = self._svariogram
            self.fold_frame_coordinate = self.fold_frame_coordinate / x_scale
            self._svariogram = None
            try:
                guess_raw = params.get("guess", None)
                if guess_raw is not None:
                    # User-supplied guess is in original space; scale w to match.
                    guess = np.array(guess_raw, dtype=float)
                    guess[-1] /= x_scale
                else:
                    guess = np.array(
                        self.initial_guess(
                            wavelength=(
                                params["wavelength"] / x_scale
                                if params.get("wavelength") is not None
                                else None
                            ),
                            reset=params.get("reset", False),
                            svariogram_parameters=params.get("svariogram_parameters", {}),
                            calculate_wavelength=params.get("calculate_wavelength", True),
                        ),
                        dtype=float,
                    )
            finally:
                self.fold_frame_coordinate = _orig_ffc
                self._svariogram = _orig_svario

            # Bounds in scaled space: w in [5%, 400%] of the scaled data range (=1.0).
            lo_w, hi_w = 0.05, 4.0
            guess[-1] = float(np.clip(guess[-1], lo_w, hi_w))
            bounds = (
                [-np.inf] * (len(guess) - 1) + [lo_w],
                [np.inf] * (len(guess) - 1) + [hi_w],
            )
            logger.info(
                f"curve_fit w bounds (scaled): [{lo_w}, {hi_w}], "
                f"initial w_scaled={guess[-1]:.4g} (~{guess[-1]*x_scale:.4g} orig)"
            )
            try:
                res = curve_fit(
                    self._function,
                    x_scaled,
                    np.tan(np.deg2rad(self.rotation_angle[mask])),
                    p0=guess,
                    bounds=bounds,
                    maxfev=5000,
                    full_output=True,
                )
                guess = res[0]
            except Exception as e:
                logger.error(f"curve_fit failed ({e}); using initial guess as fallback")

            # Scale wavelength back to original coordinate space.
            guess[-1] *= x_scale

            try:
                self.update_params(guess)
            except Exception:
                logger.error("update_params failed after fit")
                return False
            return True
        return True

    @abstractmethod
    def update_params(self, params: Union[List, npt.NDArray[np.float64]]) -> None:
        pass

    @abstractmethod
    def initial_guess(
        self,
        wavelength: Optional[float] = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict = {},
        reset: bool = False,
    ) -> np.ndarray:
        pass

    @staticmethod
    @abstractmethod
    def _function(s, *args, **kwargs):
        pass

    @property
    @abstractmethod
    def params(self) -> dict:
        pass

    def plot(self, ax=None, show_data: bool = True, **kwargs):
        if ax is None:
            import matplotlib.pyplot as plt
            _fig, ax = plt.subplots()
        if show_data and self.fold_frame_coordinate is not None:
            ax.scatter(self.fold_frame_coordinate, self.rotation_angle, c="r")
        ax.plot(self.evaluation_points, self(self.evaluation_points), **kwargs)
        return ax

    def __call__(self, s: np.ndarray) -> np.ndarray:
        return np.rad2deg(np.arctan(self._function(s, **self.params)))
