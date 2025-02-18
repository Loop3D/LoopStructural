from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile
import numpy as np
import numpy.typing as npt
from typing import Optional, Callable
from .....utils import getLogger

logger = getLogger(__name__)


class LambdaFoldRotationAngleProfile(BaseFoldRotationAngleProfile):
    def __init__(
        self,
        fn: Callable[[np.ndarray], np.ndarray],
        rotation_angle: Optional[npt.NDArray[np.float64]] = None,
        fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
    ):
        """The fold frame function using the lambda profile from Laurent 2016

        Parameters
        ----------
        rotation_angle : npt.NDArray[np.float64], optional
            the calculated fold rotation angle from observations in degrees, by default None
        fold_frame_coordinate : npt.NDArray[np.float64], optional
            fold frame coordinate scalar field value, by default None
        lambda_ : float, optional
            lambda parameter, by default 0
        """
        super().__init__(rotation_angle, fold_frame_coordinate)
        self._function = fn

    @property
    def params(self):
        return {}

    def update_params(self, params):
        pass

    def initial_guess(
        self,
        wavelength: float | None = None,
        calculate_wavelength: bool = True,
        svariogram_parameters: dict = {},
        reset: bool = False,
    ) -> np.ndarray:
        return np.array([])
