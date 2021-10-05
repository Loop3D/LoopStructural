from LoopStructural.utils import LoopImportError, LoopTypeError, LoopValueError
try:
    from LoopProjectFile import ProjectFile
except ImportError:
    raise LoopImportError("LoopProjectFile cannot be imported")

from .process_data import ProcessInputData
import numpy as np
import pandas as pd
import networkx 

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class Map2LoopProcessor(ProcessInputData):
    def __init__(self,projectife,use_thickness=None):
        if isinstance(projectife,ProjectFile) == False:
            raise LoopTypeError("projectife must be of type ProjectFile")
        self.projectife = projectife
        # super().__init__( 
        #             self.projectfile.contacts, 
        #             self.projectfile.orientations, 
        #             stratigraphic_order,
        #             thicknesses=thicknesses,
        #             fault_orientations=fault_orientations,
        #             fault_locations=fault_locations,
        #             fault_properties=fault_properties,
        #             fault_edges=list(fault_graph.edges),
        #             colours=dict(zip(groups['code'],groups['colour'])),
        #             fault_stratigraphy=None,
        #             intrusions=None,
        #             use_thickness=use_thickness,
        #             fault_edge_properties=fault_edge_properties
        #             )