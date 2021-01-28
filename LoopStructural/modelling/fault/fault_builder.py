from ..features.structural_frame_builder import StructuralFrameBuilder
class FaultBuilder(StructuralFrameBuilder):
    def __init__(self,interpolator=None,interpolators=None,**kwargs):
        """A specialised structural frame builder for building a fault

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        """
        StructuralFrameBuilder.__init__(self,interpolator,interpolators,kwargs)
    
    def create_data_from_geometry(self, 
                                data,
                                fault_center, 
                                normal_vector,
                                slip_vector,
                                influence_distance = None,
                                horizontal_radius = None,
                                vertical_radius = None):
        """Generate the required data for building a fault frame for a fault with the 
        specified parameters

        Parameters
        ----------
        data : DataFrame,
            model data
        fault_center : np.array(3)
            x,y,z coordinates of the fault center
        normal_vector : np.array(3)
            x,y,z components of normal vector to fault, single observation usually 
            average direction
        slip_vector : np.array(3)
            x,y,z components of slip vector for the fault, single observation usually 
            average direction
        influence_distance : double
            distance away from fault for the fault volume
        horizontal_radius : double
            fault extent
        vertical_radius : double
            fault volume radius in the slip direction
        """
        
        normal_vector[0] = np.sin(np.deg2rad(fault_dip_direction))
        normal_vector[1] = np.cos(np.deg2rad(fault_dip_direction))
        strike_vector = np.cross(normal_vector,slip_vector)
        if influence_distance is not None:
            fault_edges[0,:] = fault_center[:3]+normal_vector*influence_distance
            fault_edges[1,:] = fault_center[:3]-normal_vector*influence_distance
            data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = [fault_edges[0,0],fault_edges[0,1],fault_edges[0,2],f,1,0]
            data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = [fault_edges[1,0],fault_edges[1,1],fault_edges[1,2],f,-1,0]
        if horizontal_radius is not None:
            fault_tips[0,:] = fault_center[:3]+strike_vector*horizontal_radius
            fault_tips[1,:] = fault_center[:3]-strike_vector*horizontal_radius
            data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = [fault_tips[0,0],fault_tips[0,1],fault_tips[0,2],f,1,2]
            data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = [fault_tips[1,0],fault_tips[1,1],fault_tips[1,2],f,-1,2]
        if vertical_radius is not None:
            fault_depth[0,:] = fault_center[:3]+slip_vector*vertical_radius
            fault_depth[1,:] = fault_center[:3]-slip_vector*vertical_radius
            #TODO need to add data here
        
        # add strike vector to constraint fault extent

        data.loc[len(data),['X','Y','Z','feature_name','nx','ny','nz','coord']] = [fault_center[0],fault_center[1],fault_center[2],\
            f, fault_extent_vector[0], fault_extent_vector[1], fault_extent_vector[2], 2]
        self.add_data_from_data_frame(data)

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