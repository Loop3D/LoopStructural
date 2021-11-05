from LoopStructural.modelling.features import StructuralFrameBuilder

import numpy as np
import pandas as pd
from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class IntrusionBuilder(StructuralFrameBuilder):
    def __init__(self, interpolator=None, interpolators=None,
                 model=None, feature_name = None):
        """A specialised structural frame builder for building an intrusion

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        model : GeologicalModel
            reference to the model containing the fault
        """

        StructuralFrameBuilder.__init__(self,interpolator, interpolators)
        self.origin = np.array([np.nan,np.nan,np.nan])
        self.maximum = np.array([np.nan,np.nan,np.nan])
        self.model = model
        self.name = feature_name

        self.minimum_origin = self.model.bounding_box[0,:]
        self.maximum_maximum = self.model.bounding_box[1,:] 
    
    def update_geometry(self, points):
        self.origin = np.nanmin(np.array([np.min(points,axis=0),self.origin]),axis=0)
        self.maximum = np.nanmax(np.array([np.max(points,axis=0),self.maximum]),axis=0)
        self.origin[self.origin<self.minimum_origin] = self.minimum_origin[self.origin<self.minimum_origin]
        self.maximum[self.maximum>self.maximum_maximum] = self.maximum_maximum[self.maximum>self.maximum_maximum]
    
    def set_data(self, feature_data, intrusion_network_points):
        """Adds the intrusion network points as data for the intrusion frame

        Parameters
        ----------
        feature_name : string, name of the intrusion frame feature in the dataframe
        feature_data : DataFrame,
                       model data
        intrusion_network_points: numpy array [x,y,z], 
                                  outcome of IntrusionNetwork.build()
        
        """
        # Coordinate 0 - Represents growth, isovalue 0 correspond to the intrusion network surface, gradient must be provided (ix,iy,iz):
        # scaled_inet_points = self.model.scale(intrusion_network_points[:,:3])
        scaled_inet_points = intrusion_network_points[:,:3]
        coord_0_values = pd.DataFrame(scaled_inet_points, columns = ['X','Y','Z'])
        coord_0_values['val'] = 0
        coord_0_values['coord'] = 0
        coord_0_values['feature_name'] = self.name
        intrusion_frame_data = feature_data.append(coord_0_values, sort = False)
          
        self.add_data_from_data_frame(intrusion_frame_data)
        self.update_geometry(intrusion_frame_data[['X','Y','Z']].to_numpy())
        return intrusion_frame_data
        
        
    def update(self):
        for i in range(3):
            self.builders[i].update()