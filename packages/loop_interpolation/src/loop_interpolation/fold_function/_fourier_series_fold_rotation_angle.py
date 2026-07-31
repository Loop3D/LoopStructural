"""Fourier-series fold rotation-angle profile (Laurent et al., 2016)."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt
from loop_common.logging import get_logger

from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile

logger = get_logger(__name__)


class FourierSeriesFoldRotationAngleProfile(BaseFoldRotationAngleProfile):
    """Fold limb-rotation angle as a truncated Fourier series:

    ``tan(α(s)) = c0 + c1·cos(2π/w·s) + c2·sin(2π/w·s)``
    """

    def __init__(
        self,
        rotation_angle: npt.NDArray[np.float64] | None = None,
        fold_frame_coordinate: npt.NDArray[np.float64] | None = None,
        c0: float = 0.0,
        c1: float = 0.0,
        c2: float = 0.0,
        w: float = 1.0,
    ):
        super().__init__(rotation_angle, fold_frame_coordinate)
        self._c0 = c0
        self._c1 = c1
        self._c2 = c2
        self._w = w

    # --- properties with observer notification ---

    @property
    def c0(self) -> float:
        return self._c0

    @c0.setter
    def c0(self, value: float) -> None:
        self.notify_observers()
        self._c0 = value

    @property
    def c1(self) -> float:
        return self._c1

    @c1.setter
    def c1(self, value: float) -> None:
        self.notify_observers()
        self._c1 = value

    @property
    def c2(self) -> float:
        return self._c2

    @c2.setter
    def c2(self, value: float) -> None:
        self.notify_observers()
        self._c2 = value

    @property
    def w(self) -> float:
        return self._w

    @w.setter
    def w(self, value: float) -> None:
        if value <= 0:
            raise ValueError("wavelength must be > 0")
        self.notify_observers()
        self._w = value

    # --- profile function ---

    @staticmethod
    def _function(x: np.ndarray, c0: float, c1: float, c2: float, w: float) -> np.ndarray:
        return c0 + c1 * np.cos(2 * np.pi / w * x) + c2 * np.sin(2 * np.pi / w * x)

    # --- params ---

    @property
    def params(self) -> dict:
        return {"c0": self.c0, "c1": self.c1, "c2": self.c2, "w": self.w}

    @params.setter
    def params(self, params: dict) -> None:
        if "w" in params and params["w"] <= 0:
            raise ValueError("wavelength must be > 0")
        self._c0 = params["c0"]
        self._c1 = params["c1"]
        self._c2 = params["c2"]
        self._w = params["w"]
        self.notify_observers()

    def update_params(self, params: list[float] | npt.NDArray[np.float64]) -> None:
        if len(params) != 4:
            raise ValueError("params must have 4 elements: [c0, c1, c2, w]")
        self._c0 = params[0]
        self._c1 = params[1]
        self._c2 = params[2]
        self._w = params[3]
        self.notify_observers()

    def initial_guess(
        self,
        wavelength: float | None = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict | None = None,
        reset: bool = False,
    ) -> np.ndarray:
        if svariogram_parameters is None:
            svariogram_parameters = {}
        if reset:
            # Clip to ±89° before tan to avoid singularities at ±90°.
            ang = np.clip(self.rotation_angle, -89.0, 89.0)
            tan_vals = np.tan(np.deg2rad(ang))
            self.c0 = float(np.nanmean(tan_vals))
            self.c1 = 0.0
            self.c2 = float(np.nanmax(tan_vals))
            self.w = 1.0
        if calculate_wavelength:
            self.w = self.estimate_wavelength(svariogram_parameters=svariogram_parameters)
        if wavelength is not None:
            self.w = wavelength
        return np.array([self.c0, self.c1, self.c2, self.w])

    def calculate_misfit(
        self, s: np.ndarray, rotation_angle: np.ndarray
    ) -> np.ndarray:
        return np.tan(np.deg2rad(rotation_angle)) - np.tan(np.deg2rad(self.__call__(s)))
