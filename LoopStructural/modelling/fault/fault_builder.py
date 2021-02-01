from ..features.structural_frame_builder import StructuralFrameBuilder

import numpy as np
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

        StructuralFrameBuilder.__init__(self,interpolator,interpolators,**kwargs)
        self.origin = np.array([np.nan,np.nan,np.nan])
        self.maximum = np.array([np.nan,np.nan,np.nan])

    def update_geometry(self,points):
        self.origin = np.nanmin(np.array([np.min(points,axis=0),self.origin]),axis=0)
        self.maximum = np.nanmax(np.array([np.max(points,axis=0),self.maximum]),axis=0)

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
        strike_vector = np.cross(normal_vector,slip_vector)
        fault_edges = np.zeros((2,3))
        fault_tips = np.zeros((2,3))
        fault_depth = np.zeros((2,3))
        if fault_center is not None:
            if influence_distance is not None:
                fault_edges[0,:] = fault_center[:3]+normal_vector*influence_distance
                fault_edges[1,:] = fault_center[:3]-normal_vector*influence_distance
                self.update_geometry(fault_edges)
                data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                    [fault_edges[0,0],fault_edges[0,1],fault_edges[0,2],self.name,1,0]
                data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                    [fault_edges[1,0],fault_edges[1,1],fault_edges[1,2],self.name,-1,0]
            if horizontal_radius is not None:
                fault_tips[0,:] = fault_center[:3]+strike_vector*horizontal_radius
                fault_tips[1,:] = fault_center[:3]-strike_vector*horizontal_radius
                self.update_geometry(fault_tips)
                data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                    [fault_tips[0,0],fault_tips[0,1],fault_tips[0,2],self.name,1,2]
                data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                    [fault_tips[1,0],fault_tips[1,1],fault_tips[1,2],self.name,-1,2]
            if vertical_radius is not None:
                fault_depth[0,:] = fault_center[:3]+slip_vector*vertical_radius
                fault_depth[1,:] = fault_center[:3]-slip_vector*vertical_radius
                self.update_geometry(fault_depth)
                #TODO need to add data here
            data.loc[len(data),['X','Y','Z','feature_name','nx','ny','nz','coord']] =\
                [fault_center[0],fault_center[1],fault_center[2],self.name,slip_vector[0],slip_vector[1],slip_vector[2],1]
        # add strike vector to constraint fault extent
            data.loc[len(data),['X','Y','Z','feature_name','nx','ny','nz','coord']] = [fault_center[0],fault_center[1],fault_center[2],\
                self.name, strike_vector[0], strike_vector[1], strike_vector[2], 2]
        self.add_data_from_data_frame(data)
    def set_mesh_geometry(self,buffer):
        """set the mesh geometry

        Parameters
        ----------
        buffer : double 
            percentage of length to add to edges
        """
        length = self.maximum-self.origin

        # for builder in self.builders:
        # all three coordinates share the same support
        self.builders[0].set_interpolation_geometry(self.origin-length*buffer,self.maximum-length*buffer)
            
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