from ..features.structural_frame_builder import StructuralFrameBuilder
class FaultBuilder(StructuralFrameBuilder):
    def __init__(self,interpolator=None,interpolators=None,**kwargs):
        StructuralFrameBuilder.__init__(self,interpolator,interpolators,kwargs)
    
    def create_data_from_geometry(self, 
                                fault_center, 
                                fault_dip_direction,
                                influence_distance,
                                horizontal_radius,
                                vertical_radius):
        
        normal_vector[0] = np.sin(np.deg2rad(fault_dip_direction))
        normal_vector[1] = np.cos(np.deg2rad(fault_dip_direction))
        strike_vector[0] = normal_vector[1]
        strike_vector[1] = -normal_vector[0]
        slip_vector[2]=1
        fault_edges[0,:] = fault_center[:3]+normal_vector*influence_distance
        fault_edges[1,:] = fault_center[:3]-normal_vector*influence_distance
        fault_tips[0,:] = fault_center[:3]+strike_vector*horizontal_radius
        fault_tips[1,:] = fault_center[:3]-strike_vector*horizontal_radius
        fault_depth[0,:] = fault_center[:3]+slip_vector*vertical_radius
        fault_depth[1,:] = fault_center[:3]-slip_vector*vertical_radius
        fault_locations.loc[len(fault_locations),['X','Y','Z','feature_name','val','coord']] = [fault_edges[0,0],fault_edges[0,1],fault_edges[0,2],f,1,0]
        fault_locations.loc[len(fault_locations),['X','Y','Z','feature_name','val','coord']] = [fault_edges[1,0],fault_edges[1,1],fault_edges[1,2],f,-1,0]
        fault_locations.loc[len(fault_locations),['X','Y','Z','feature_name','val','coord']] = [fault_tips[0,0],fault_tips[0,1],fault_tips[0,2],f,1,2]
        fault_locations.loc[len(fault_locations),['X','Y','Z','feature_name','val','coord']] = [fault_tips[1,0],fault_tips[1,1],fault_tips[1,2],f,-1,2]
        # add strike vector to constraint fault extent
        fault_orientations.loc[len(fault_orientations),['X','Y','Z','feature_name','nx','ny','nz','coord']] = [fault_center[0],fault_center[1],fault_center[2],f, fault_centers[3]-90,2]
    
    def add_splay(self,splayregion,splay):
        for i in range(3):
            # work out the values of the nodes where we want hard
            # constraints
            idc = np.arange(0, interpolator.support.n_nodes)[
                kwargs['splayregion'](interpolator.support.nodes)]
            val = kwargs['splay'][i].evaluate_value(
                interpolator.support.nodes[
                kwargs['splayregion'](interpolator.support.nodes), :])
            mask = ~np.isnan(val)
            fault_frame_builder[i].interpolator.add_equality_constraints(
                idc[mask], val[mask])