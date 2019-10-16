from LoopStructural.modelling.core.geological_points import GPoint, IPoint
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator, GeologicalFeature
from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.modelling.features.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fault.fault_segment import FaultSegment

import numpy as np
import networkx as nx
import pandas as pd

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
        # self.bounding_box[1, :] = self.maximum-self.origin
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
        # self.data['X'] -= self.origin[0]
        # self.data['Y'] -= self.origin[1]
        # self.data['Z'] -= self.origin[2]
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

    def get_interpolator(self, interpolatortype = 'PLI', nelements = 5e2, buffer=1e-1,**kwargs):
        """
        Returns an interpolator given the arguments
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
        interpolator = None
        bb = np.copy(self.bounding_box)
        bb[0,:] -= buffer
        bb[1,:] += buffer
        if interpolatortype == "PLI":
            mesh = TetMesh()

            mesh.setup_mesh(bb, n_tetra=nelements,)
            return PLI(mesh)

        if interpolatortype == 'FDI':
            # number of elements should be divided roughly to match the shape
            ratio = bb[1,:] / np.sum(bb[1,:])
            ratio *= nelements
            ratio = ratio.astype(int)
            step_vector = 1. / np.max(ratio)
            grid = StructuredGrid(nsteps=ratio, step_vector=step_vector)
            return FDI(grid)
        logger.warning("No interpolator")
        return interpolator

    def create_and_add_conformable_series(self, series_surface_data, **kwargs):
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
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        #
        fault_frame_builder = StructuralFrameBuilder(interpolator,**kwargs)
        # add data
        # if there is no fault slip data then
        fault_frame = fault_frame_builder.build(**kwargs)
        #
        for f in reversed(self.features):
            if f.type is 'unconformity':
                fault_frame[:].add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break
        fault = FaultSegment(fault_frame,**kwargs)
        self.features.append(fault)
        return fault

    def add_fold(self, fold_data, **kwargs):
        self.graph.add_node(fold, name=fold.name)

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
