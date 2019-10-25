from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.modelling.features.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fault.fault_segment import FaultSegment

import numpy as np
import networkx as nx
import pandas as pd
from scipy.optimize import minimize
import logging
logger = logging.getLogger(__name__)


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes the rescaling
    of the model between 0 and 1.
    """
    def __init__(self, origin, maximum):
        """

        Parameters
        ----------
        origin - numpy array specifying the origin of the model
        maximum - numpy array specifying the maximum extent of the model
        """
        self.graph = nx.DiGraph()
        self.features = []
        self.data = None


        # we want to rescale the model area so that the maximum length is
        # 1
        self.origin = np.array(origin)

        self.maximum = np.array(maximum)
        lengths = self.maximum - self.origin
        # self.scale_factor = np.max(lengths)

        self.bounding_box = np.zeros((2, 3))
        self.bounding_box[1, :] = self.maximum-self.origin
        # self.bounding_box /= self.scale_factor
        self.bounding_box[0,:] = self.origin
        self.bounding_box[1,:] = self.maximum
    def set_model_data(self, data):
        """
        Set the data array for the model
        Parameters
        ----------
        data - pandas data frame with column headers corresponding to the
         type, X, Y, Z, nx, ny, nz, val, strike, dip, dip_dir, plunge, plunge_dir, azimuth

        Returns
        -------

        Note
        ----
        Type can be any unique identifier for the feature the data point 'eg' 'S0', 'S2', 'F1_axis'
        it is then used by the create functions to get the correct data
        """
        self.data = data
        self.data['X'] -= self.origin[0]
        self.data['Y'] -= self.origin[1]
        self.data['Z'] -= self.origin[2]
        # self.data['X'] /= self.scale_factor
        # self.data['Y'] /= self.scale_factor
        # self.data['Z'] /= self.scale_factor

    def extend_model_data(self, newdata):
        """
        Extends the data frame
        Parameters
        ----------
        newdata - pandas data frame

        Returns
        -------

        """
        self.data.append(newdata)

    def get_interpolator(self, interpolatortype = 'PLI', nelements = 5e2, buffer=0.2,**kwargs):
        """
        Returns an interpolator given the arguments, also constructs a support for a discrete
        interpolator
        Parameters
        ----------
        interpolatortype - string
            define the interpolator type
        nelements - int
            number of elements in the interpolator
        buffer - double or numpy array 3x1
            value(s) between 0,1 specifying the buffer around the bounding box

        kwargs - no kwargs used, this just catches any additional arguments

        Returns
        -------

        """
        # get an interpolator for 
        interpolator = None
        bb = np.copy(self.bounding_box)
        bb[0,:] -= buffer*(bb[1,:]-bb[0,:])
        bb[1,:] += buffer*(bb[1,:]-bb[0,:])
        if interpolatortype == "PLI":
            mesh = TetMesh()

            mesh.setup_mesh(bb, n_tetra=nelements,)
            return PLI(mesh)

        if interpolatortype == 'FDI':
            # find the volume of one element
            ele_vol = bb[1,0]*bb[1,1]*bb[1,2] / nelements
            # calculate the relative lengths of the volume (x+y+z = 1)
            step_vector = bb[1,:] / np.sum(bb[1,:])
            # length of element ratio*scale = element_volume / cuberoot(l1*l2*l3)
            scale = (ele_vol / (step_vector[0]*step_vector[1]*step_vector[2]))**(1/3)
            step_vector*=scale

            # ratio = ratio.astype(int)
            # round up nsteps = length of volume / cell size
            nsteps = np.ceil(bb[1,0]/step_vector).astype(int)
            grid = StructuredGrid(origin=bb[0,:],nsteps=nsteps, step_vector=step_vector)
            return FDI(grid)
        logger.warning("No interpolator")
        return interpolator

    def create_and_add_conformable_foliation(self, series_surface_data, **kwargs):
        """
        Parameters
        ----------
        series_surface_data - string corresponding to the type in the data frame column
        kwargs

        Returns
        -------

        """
        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureInterpolator(interpolator,name=series_surface_data)
        # add data
        series_data = self.data[self.data['type'] == series_surface_data]
        series_builder.add_data_from_data_frame(series_data)
        for f in reversed(self.features):
            if f.type == 'fault':
                series_builder.add_fault(f)
            if f.type == 'unconformity':
                break
        # build feature
        series_feature = series_builder.build(**kwargs)
        series_feature.type = 'series'
        # see if any unconformities are above this feature if so add region
        for f in reversed(self.features):
            if f.type is 'unconformity':
                series_feature.add_region(lambda pos: f.evaluate_value(pos) < 0)
                break
        self.features.append(series_feature)
        return series_feature

    def create_and_add_unconformity(self, unconformity_surface_data,**kwargs):
        """

        Parameters
        ----------
        unconformity_surface_data

        Returns
        -------

        """
        interpolator = self.get_interpolator(**kwargs)
        unconformity_feature_builder = GeologicalFeatureInterpolator(interpolator, name=unconformity_surface_data)
        # add data
        unconformity_data = self.data[self.data['type'] == unconformity_surface_data]

        unconformity_feature_builder.add_data_from_data_frame(unconformity_data)
        # look through existing features if there is a fault before an unconformity
        # then add to the feature, once we get to an unconformity stop
        for f in reversed(self.features):
            if f.type == 'fault':
                unconformity_feature_builder.add_fault(f)
            if f.type == 'unconformity':
                break

        # build feature
        uc_feature = unconformity_feature_builder.build(**kwargs)
        uc_feature.type = 'unconformity'

        # iterate over existing features and add the unconformity as a region so the feature is only
        # evaluated where the unconformity is positive
        for f in self.features:
            f.add_region(lambda pos: uc_feature.evaluate_value(pos) >= 0)

        # see if any unconformities are above this feature if so add region
        for f in reversed(self.features):
            if f.type == 'unconformity':
                uc_feature.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break

        self.features.append(uc_feature)
        return uc_feature

    def create_and_add_fault(self, fault_surface_data, displacement, **kwargs):
        """

        Parameters
        ----------
        fault_surface_data - string
            name of the fault surface data in the dataframe
        displacement - displacement magnitude
        kwargs - additional kwargs for Fault and interpolators

        Returns
        -------

        """
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        #
        fault_frame_builder = StructuralFrameBuilder(interpolator,name=fault_surface_data,**kwargs)
        # add data
        fault_frame_data = self.data[self.data['type'] == fault_surface_data]
        fault_frame_builder.add_data_from_data_frame(fault_frame_data)
        # if there is no fault slip data then we could find the strike of the fault and build
        # the second coordinate
        fault_frame = fault_frame_builder.build(**kwargs)
        # if we add a region to the fault then the fault operator doesn't work but for visualisation
        # we want to add a region!
        # for f in reversed(self.features):
        #     if f.type is 'unconformity':
        #         fault_frame[0].add_region(lambda pos: f.evaluate_value(pos) <= 0)
        #         fault_frame[1].add_region(lambda pos: f.evaluate_value(pos) <= 0)
        #         fault_frame[2].add_region(lambda pos: f.evaluate_value(pos) <= 0)
        #         break
        fault = FaultSegment(fault_frame,displacement=displacement,**kwargs)
        self.features.append(fault)
        return fault

    def create_and_add_folded_foliation(self, folded_foliation_data, fold_frame, **kwargs):

        pass

    def add_feature(self, feature, name):
        self.features[name] = feature

    def rescale(self, points):
        points*=self.scale_factor
        points+=self.origin
        return points

    def view(self):

        pass
    def voxet(self, nsteps = (50, 50, 25)):
        return {'bounding_box': self.bounding_box, 'nsteps': nsteps}
