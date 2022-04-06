from enum import IntEnum
from LoopStructural.utils import getLogger


class FeatureType(IntEnum):
    """
    Enum for the different interpolator types

    1-9 should cover interpolators with supports
    9+ are data supported
    """

    BASE = 0
    INTERPOLATED = 1
    STRUCTURALFRAME = 2
    REGION = 3
    FOLDEDFEATURE = 4
    ANALYTICALFEATURE = 5
    LAMBDAFEATURE = 6


from ._base_geological_feature import BaseFeature
from .geological_feature import GeologicalFeature
from .lambda_geological_feature import LambdaGeologicalFeature
from .geological_feature_builder import GeologicalFeatureInterpolator
from .region_feature import RegionFeature
from .structural_frame import StructuralFrame
from .structural_frame_builder import StructuralFrameBuilder
from .unconformity_feature import UnconformityFeature
from .analytical_feature import AnalyticalGeologicalFeature
