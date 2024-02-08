from enum import IntEnum


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
    DOMAINFAULT = 10
    INACTIVEFAULT = 11



# from .builders._geological_feature_builder import GeologicalFeatureBuilder

