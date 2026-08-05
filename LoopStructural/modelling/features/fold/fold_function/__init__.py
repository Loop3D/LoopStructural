from __future__ import annotations

from enum import Enum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._fourier_series_fold_rotation_angle import FourierSeriesFoldRotationAngleProfile
from ._lambda_fold_rotation_angle import LambdaFoldRotationAngleProfile
from ._trigo_fold_rotation_angle import TrigoFoldRotationAngleProfile


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
    rotation_angle: npt.NDArray[np.float64] | None = None,
    fold_frame_coordinate: npt.NDArray[np.float64] | None = None,
    **kwargs,
):
    return fold_rotation_type.value(rotation_angle, fold_frame_coordinate, **kwargs)
