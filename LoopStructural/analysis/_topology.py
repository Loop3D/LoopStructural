import numpy as np
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
    # topology_matrix[:] = -9999
    for i, f in enumerate(model.faults):
        topology_matrix[:,i] = f.evaluate(xyz)
    topology_matrix[topology_matrix==0] = -1
    topology_matrix[np.isnan(topology_matrix)] =0 
    # topology_matrix[topology_matrix==0] = -1
    # topology_matrix[topology_matrix==-9999] = 0
    return topology_matrix.astype(int)