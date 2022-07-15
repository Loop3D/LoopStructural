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
    FOLDED = 4
    ANALYTICAL = 5
    LAMBDA = 6
    UNCONFORMITY = 7
    INTRUSION = 8
    FAULT = 9


from ._base_geological_feature import BaseFeature
from ._geological_feature import GeologicalFeature
from ._lambda_geological_feature import LambdaGeologicalFeature

# from .builders._geological_feature_builder import GeologicalFeatureBuilder
from ._structural_frame import StructuralFrame
from ._cross_product_geological_feature import CrossProductGeologicalFeature

from ._unconformity_feature import UnconformityFeature
from ._analytical_feature import AnalyticalGeologicalFeature
from ._structural_frame import StructuralFrame
