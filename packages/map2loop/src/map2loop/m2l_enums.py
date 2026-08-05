from enum import IntEnum


class Datatype(IntEnum):
    GEOLOGY = 0
    STRUCTURE = 1
    FAULT = 2
    FOLD = 3
    DTM = 4
    FAULT_ORIENTATION = 5


class Datastate(IntEnum):
    UNNAMED = 0
    UNLOADED = 1
    LOADED = 2
    REPROJECTED = 3
    CLIPPED = 4
    COMPLETE = 5
    ERRORED = 9


class ErrorState(IntEnum):
    NONE = 0
    URLERROR = 1
    CONFIGERROR = 2


class VerboseLevel(IntEnum):
    NONE = 0
    TEXTONLY = 1
    ALL = 2
