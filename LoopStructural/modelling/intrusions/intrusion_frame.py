import numpy as np

from LoopStructural.modelling.features.structural_frame import StructuralFrame

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class IntrusionFrame(StructuralFrame):
    def __init__(self, name, features): 
        StructuralFrame.__init__(self, name, features)
        logger.warning("There is an error")

        

