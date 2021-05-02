def _calculate_fault_topology_matrix(self, model, xyz):
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
    topology_matrix = np.zeros((xyz.shape,len(model.faults)),dtype=int)
    for i, f in enumerate(model.faults):
        topology_matrix[:,i].evaluate(xyz)