from ._trigo_fold_rotation_angle import TrigoFoldRotationAngleProfile
from ._fourier_series_fold_rotation_angle import FourierSeriesFoldRotationAngleProfile
from enum import Enum


class FoldRotationType(Enum):
    TRIGONOMETRIC = TrigoFoldRotationAngleProfile
    FOURIER_SERIES = FourierSeriesFoldRotationAngleProfile
    # ADDITIONAL = AdditionalFoldRotationAngle

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name
