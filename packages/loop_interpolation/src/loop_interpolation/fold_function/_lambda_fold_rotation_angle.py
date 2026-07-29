"""Lambda (arbitrary callable) fold rotation-angle profile."""

from typing import Callable, Optional

import numpy as np
import numpy.typing as npt

from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile


class LambdaFoldRotationAngleProfile(BaseFoldRotationAngleProfile):
    """Fold rotation-angle profile backed by an arbitrary callable.

    ``__call__(s)`` simply delegates to the supplied function and returns
    degrees (the function is expected to return degrees directly).
    """

    def __init__(
        self,
        fn: Callable[[np.ndarray], np.ndarray],
        rotation_angle: Optional[npt.NDArray[np.float64]] = None,
        fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
    ):
        super().__init__(rotation_angle, fold_frame_coordinate)
        self._fn = fn

    # Override __call__ — the base computes arctan(_function(...)), but for a
    # lambda profile the function already returns angles in degrees.
    def __call__(self, s: np.ndarray) -> np.ndarray:
        return self._fn(s)

    @staticmethod
    def _function(s, *args, **kwargs):
        raise NotImplementedError("LambdaFoldRotationAngleProfile uses a supplied callable")

    @property
    def params(self) -> dict:
        return {}

    def update_params(self, params) -> None:
        pass

    def initial_guess(
        self,
        wavelength: Optional[float] = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict = {},
        reset: bool = False,
    ) -> np.ndarray:
        return np.array([])
