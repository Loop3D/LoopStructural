from ._trigo_fold_rotation_angle import TrigoFoldRotationAngleProfile
from ._fourier_series_fold_rotation_angle import FourierSeriesFoldRotationAngleProfile
from enum import Enum
from typing import Optional
import numpy.typing as npt
import numpy as np


class FoldRotationType(Enum):
    TRIGONOMETRIC = TrigoFoldRotationAngleProfile
    FOURIER_SERIES = FourierSeriesFoldRotationAngleProfile
    # ADDITIONAL = AdditionalFoldRotationAngle

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


def get_fold_rotation_profile(
    fold_rotation_type,
    rotation_angle: Optional[npt.NDArray[np.float64]] = None,
    fold_frame_coordinate: Optional[npt.NDArray[np.float64]] = None,
    **kwargs,
):
    return fold_rotation_type.value(rotation_angle, fold_frame_coordinate, **kwargs)
