class LoopStructuralConfig:
    """
    Class to store configuration settings.
    """

    __splay_fault_threshold = 30
    __experimental = False
    __default_interpolator = "FDI"
    __default_nelements = 1e4
    __default_solver = "cg"

    # @property
    # def experimental():
    #     return __experimental

    # @experimental.setter
    # def experimental(self, value):
    #     __experimental = value
