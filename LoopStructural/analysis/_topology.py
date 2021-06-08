import numpy as np
from LoopStructural.utils import getLogger
logger = getLogger(__name__)

def calculate_fault_topology_matrix(model, xyz=None):
    """Calculate fault ellipsoid and hw/fw

    Parameters
    ----------
    model : GeologicalModel
        the model containing the faults
    xyz : np.array
        xyz locations in model coordinates
    
    Returns
    -------
    topology_matrix : np.array
        matrix containing nan (outside), 0 (footwall), 1 (hangingwall)
    """
    if xyz is None:
        xyz = model.regular_grid(rescale=False,shuffle=False)
    topology_matrix = np.zeros((xyz.shape[0],len(model.faults)))
    topology_matrix[:] = np.nan
    for i, f in enumerate(model.faults):
        topology_matrix[f.inside_volume(xyz),i] = f.evaluate(xyz[f.inside_volume(xyz),:])
    topology_matrix[topology_matrix==0] = -1
    topology_matrix[np.isnan(topology_matrix)] =0 
    return topology_matrix.astype(int)


    