from ..features.structural_frame_builder import StructuralFrameBuilder

import numpy as np
class FaultBuilder(StructuralFrameBuilder):
    def __init__(self,interpolator=None,interpolators=None,model=None,fault_bounding_box_buffer=0.2,**kwargs):
        """A specialised structural frame builder for building a fault

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        model : GeologicalModel
            reference to the model containing the fault
        fault_bounding_box_buffer: float, default 0.2
            the maximum area around the model domain that a fault is modelled. For high displacement faults this
            may need to be large, smaller values will be result in fewer degrees of freedom = quicker interpolation
        """

        StructuralFrameBuilder.__init__(self,interpolator,interpolators,**kwargs)
        self.origin = np.array([np.nan,np.nan,np.nan])
        self.maximum = np.array([np.nan,np.nan,np.nan])
        self.model = model
        # define a maximum area to mesh adding buffer to model
        buffer = .2
        self.minimum_origin = self.model.bounding_box[0,:] - buffer*(self.model.bounding_box[1,:]-self.model.bounding_box[0,:])
        self.maximum_maximum = self.model.bounding_box[1,:] + buffer*(self.model.bounding_box[1,:]-self.model.bounding_box[0,:])

    def update_geometry(self,points):
        
        self.origin = np.nanmin(np.array([np.min(points,axis=0),self.origin]),axis=0)
        self.maximum = np.nanmax(np.array([np.max(points,axis=0),self.maximum]),axis=0)
        self.origin[self.origin<self.minimum_origin] = self.minimum_origin[self.origin<self.minimum_origin]
        self.maximum[self.maximum>self.maximum_maximum] = self.maximum_maximum[self.maximum>self.maximum_maximum]

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
        normal_vector/=np.linalg.norm(normal_vector)
        slip_vector/=np.linalg.norm(slip_vector)
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
                    [fault_center[0],fault_center[1],fault_center[2],self.name,0,2]
                # data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                #     [fault_tips[1,0],fault_tips[1,1],fault_tips[1,2],self.name,-1,2]
                strike_vector /= horizontal_radius
            if vertical_radius is not None:
                fault_depth[0,:] = fault_center[:3]+slip_vector*vertical_radius
                fault_depth[1,:] = fault_center[:3]-slip_vector*vertical_radius
                data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                    [fault_center[0],fault_center[1],fault_center[2],self.name,0,1]
                # data.loc[len(data),['X','Y','Z','feature_name','val','coord']] = \
                #     [fault_depth[1,0],fault_depth[1,1],fault_depth[1,2],self.name,-1,1]
                self.update_geometry(fault_depth)
                #TODO need to add data here
                slip_vector /= vertical_radius

            data.loc[len(data),['X','Y','Z','feature_name','nx','ny','nz','val','coord']] =\
                [fault_center[0],fault_center[1],fault_center[2],self.name,slip_vector[0],slip_vector[1],slip_vector[2],0,1]
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
        # origin = self.builders[0].interpolator.support.origin
        # maximum = self.builders[0].interpolator.support.maximum#set_interpolation_geometry
        # if origin[2]>self.origin[2]:
        #     origin[2]=self.origin[2]
        # if maximum[2]<self.maximum[2]:
        #     maximum[2]=self.maximum[2]
        # self.builders[0].set_interpolation_geometry(origin,maximum)
        # for builder in self.builders:
        # all three coordinates share the same support
        self.builders[0].set_interpolation_geometry(self.origin-length*buffer,self.maximum+length*buffer)
            
    def add_splay(self,splayregion,splay):
        # for i in range(3):
        #     # work out the values of the nodes where we want hard
        #     # constraints
        #     idc = np.arange(0, interpolator.support.n_nodes)[
        #         kwargs['splayregion'](interpolator.support.nodes)]
        #     val = kwargs['splay'][i].evaluate_value(
        #         interpolator.support.nodes[
        #         kwargs['splayregion'](interpolator.support.nodes), :])
        #     mask = ~np.isnan(val)
        #     fault_frame_builder[i].interpolator.add_equality_constraints(
        #         idc[mask], val[mask])
        pass
    def update(self):
        for i in range(3):
            self.builders[i].update()