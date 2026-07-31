"""Fold rotation-angle profile implementations."""
from __future__ import annotations

from enum import Enum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._base_fold_rotation_angle import BaseFoldRotationAngleProfile
from ._fourier_series_fold_rotation_angle import FourierSeriesFoldRotationAngleProfile
from ._lambda_fold_rotation_angle import LambdaFoldRotationAngleProfile

__all__ = [
    "BaseFoldRotationAngleProfile",
    "FoldRotationType",
    "FourierSeriesFoldRotationAngleProfile",
    "LambdaFoldRotationAngleProfile",
    "get_fold_rotation_profile",
]


class FoldRotationType(Enum):
    FOURIER_SERIES = FourierSeriesFoldRotationAngleProfile

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name


def get_fold_rotation_profile(
    fold_rotation_type: FoldRotationType,
    rotation_angle: npt.NDArray[np.float64] | None = None,
    fold_frame_coordinate: npt.NDArray[np.float64] | None = None,
    **kwargs,
) -> BaseFoldRotationAngleProfile:
    return fold_rotation_type.value(rotation_angle, fold_frame_coordinate, **kwargs)
