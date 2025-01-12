from enum import IntEnum


class FeatureType(IntEnum):
    """ """

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
    DOMAINFAULT = 10
    INACTIVEFAULT = 11
    ONLAPUNCONFORMITY = 12


# from .builders._geological_feature_builder import GeologicalFeatureBuilder
from ._base_geological_feature import BaseFeature
from ._geological_feature import GeologicalFeature
from ._lambda_geological_feature import LambdaGeologicalFeature

# from .builders._geological_feature_builder import GeologicalFeatureBuilder
from ._structural_frame import StructuralFrame
from ._cross_product_geological_feature import CrossProductGeologicalFeature

from ._unconformity_feature import UnconformityFeature
from ._analytical_feature import AnalyticalGeologicalFeature
from ._projected_vector_feature import ProjectedVectorFeature
