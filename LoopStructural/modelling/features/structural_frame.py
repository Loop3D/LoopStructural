from LoopStructural.modelling.core.geological_points import GPoint, IPoint, TPoint
from LoopStructural.modelling.features.geological_feature import GeologicalFeature, GeologicalFeatureInterpolator
from LoopStructural.modelling.features.cross_product_geological_feature import CrossProductGeologicalFeature
from LoopStructural.supports.scalar_field import ScalarField

import numpy as np

import logging
logger = logging.getLogger(__name__)

class StructuralFrame:
    def __init__(self, name, features):
        """
        Structural frame is a curvilinear coordinate system defined by structural
        observations associated with a fault or fold.

        Parameters
        ----------
        name - name of the structural frame
        features - list of features to build the frame with
        """
        self.name = name
        self.features = features
        self.data = None

    def __getitem__(self, item):
        """

        Parameters
        ----------
        item index of feature to access

        Returns
        -------
        the structural frame geological feature
        """
        return self.features[item]

    def get_feature(self, i):
        """
        Return the ith feature

        Parameters
        ----------
        i

        Returns
        -------

        """
        return self.features[i]

    def set_data(self, data):
        """
        Associate data with structural frame

        Parameters
        ----------
        data

        Returns
        -------

        """
        self.data = data

    def evaluate_value(self, evaluation_points, i=None):
        """
        Evaluate the value of the structural frame for the points.
        Can optionally only evaluate one coordinate

        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        if i is not None:
            self.features[i].support.evaluate_value(evaluation_points)
        return (self.features[0].support.evaluate_value(evaluation_points),
                self.features[1].support.evaluate_value(evaluation_points),
                self.features[2].support.evaluate_value(evaluation_points))

    def evaluate_gradient(self, evaluation_points, i=None):
        """
        Evaluate the gradient of the structural frame.
        Can optionally only evaluate the ith coordinate

        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        if i is not None:
            return self.features[i].support.evaluate_gradient(evaluation_points)
        return (self.features[0].support.evaluate_gradient(evaluation_points),
                self.features[1].support.evaluate_gradient(evaluation_points),
                self.features[2].support.evaluate_gradient(evaluation_points))


class StructuralFrameBuilder:
    def __init__(self, interpolator = None, interpolators = None, **kwargs):
        """
        Class for building a structural frame - has functions to set up the interpolator with
        data and also orthogonality constraints. Can build a generic structural frame or a
        subclass of the structural frame if the kwarg `frame` is specified

        Parameters
        ----------
        interpolator - a template interpolator for the frame
        kwargs
        """
        self.support = None
        self.fault_event = None
        self.name = 'Undefined'
        #self.region = 'everywhere'
        self.builders = []
        # if 'support' in kwargs:
        #     self.support = kwargs['support']
        #     del kwargs['support']
        # if 'mesh' in kwargs:
        #     self.support = kwargs['mesh']
        #     del kwargs['mesh']
        if 'name' in kwargs:
            self.name = kwargs['name']
            kwargs.pop('name')
        # if 'region' in kwargs:
        #     self.region = kwargs['region']
        # if self.support is None:
        #     self.support = interpolator.support
        self.data = [ [] , [] ,[]]
        # list of interpolators
        # self.interpolators = []
        # Create the interpolation objects by copying the template
        
        if interpolators is None:
            interpolators = []
            interpolators.append(interpolator)
            interpolators.append(interpolator.copy())
            interpolators.append(interpolator.copy())
        if interpolators is None and interpolator is None:
            raise BaseException
        # self.builders
        self.builders.append(GeologicalFeatureInterpolator(interpolators[0], name = self.name + '_0', **kwargs))#,region=self.region))
        self.builders.append(GeologicalFeatureInterpolator(interpolators[1], name = self.name + '_1', **kwargs))#,region=self.region))
        self.builders.append(GeologicalFeatureInterpolator(interpolators[2], name = self.name + '_2', **kwargs))#,region=self.region))

    def __getitem__(self, item):
        return self.builders[item]

    def add_data_from_data_frame(self, data_frame):
        """
        extract the data for a fault from a data frame

        Parameters
        ----------
        data_frame

        Returns
        -------

        """
        if 'X' not in data_frame.columns or 'X' not in data_frame.columns or 'X' not in data_frame.columns:
            logger.error("No location in data frame")
            return
        for i, r in data_frame.iterrows():

            if np.isnan(r['X']) or np.isnan(r['Y']) or np.isnan(r['Z']):
                logger.debug("X,Y,Z all NAN")
                continue
            pos = r[['X', 'Y', 'Z']]
            # by default add everything to the first coordinate of the structural frame
            coord = 0
            if 'coord' in data_frame.columns and ~np.isnan(r['coord']):
                coord = int(r['coord'])
            if 'val' in data_frame.columns and ~np.isnan(r['val']):
                self.add_point(pos, r['val'],coord=coord)
            if 'strike' in data_frame.columns and 'dip' in data_frame.columns and  \
                    ~np.isnan(r['strike']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(r['polarity']):
                    polarity = r['polarity']
                self.add_strike_and_dip(pos, r['strike'], r['dip'], polarity=polarity,
                                        coord=coord)
            if 'azimuth' in data_frame.columns and 'dip' in data_frame.columns and \
                ~np.isnan(r['azimuth']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(r['polarity']):
                    polarity = r['polarity']
                self.add_plunge_and_plunge_dir(pos,r['dip'],r['azimuth'],
                                               polarity=polarity,coord=coord)
            if 'nx' in data_frame.columns and 'ny' in data_frame.columns and 'nz' in data_frame.columns and \
                    ~np.isnan(r['nx']) and ~np.isnan(r['ny'])and ~np.isnan(r['nz']):
                 self.add_planar_constraint(pos,r[['nx','ny','nz']],
                         coord=coord)

    def add_strike_dip_and_value(self, pos, strike, dip, val, polarity = 1, coord = None, itype= None):
        """
        Add a planar measurement and value to the interpolator for a coordinate of the fold frame

        Parameters
        ----------
        pos: numpy array
            position of measurement
        strike: double
            strike of plane
        dip: double
            dip of plane
        val:
            scalar value of plane
        polarity: int
            -1 is inverse polarity default 1 is standard
        coord: int
            index of fold frame
        itype:
            depreciated

        Returns
        -------

        """
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.builders[coord].add_strike_and_dip(pos, strike, dip, polarity)
        self.builders[coord].add_point(pos,val)

    def add_point(self, pos, val,  coord = None, itype= None):
        """

        Parameters
        ----------
        pos: numpy array
            position of measurement
        val: double
            scalar value
        coord: int
            structural frame coordinate
        itype:
            depreciated

        Returns
        -------

        """
        if itype == 'gx':
            logger.warning("itype will be removed, use coord instead")
            coord = 0
        if itype == 'gy':
            logger.warning("itype will be removed, use coord instead")
            coord = 1
        if itype == 'gz':
            logger.warning("itype will be removed, use coord instead")
            coord = 2
        self.builders[coord].add_point(pos,val)

    def add_planar_constraint(self, pos, val,  coord = None, itype= None):
        """

        Parameters
        ----------
        pos: numpy array
            location of observations (x,y,z)
        val: numpy array
            normal vector to plane [nx,ny,nz]
        coord: int
            index of structural frame
        itype:
            DO NOT USE depreciated

        Returns
        -------

        """
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2

        self.builders[coord].add_planar_constraint(pos, val)

    def add_plunge_and_plunge_dir(self, pos, plunge, plunge_dir, polarity =1 , coord = None, itype= None):
        """

        Parameters
        ----------
        pos: numpy array
            location of observation [x,y,z]
        plunge: double
            plunge
        plunge_dir: double
            plunge direction
        polarity: int
            polarity of vector
        coord: int
            coordinate of structural frame
        itype:
            depreciated

        Returns
        -------

        """
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.builders[coord].add_plunge_and_plunge_dir(pos,plunge,plunge_dir, polarity)

    def add_strike_and_dip(self, pos, strike, dip, polarity = 1, coord = None, itype= None):
        """

        Parameters
        ----------
        pos: numpy array
            position of observations [x,y,z]
        strike: double
            strike of plane
        dip: double
            dip of plane
        polarity: int
            polarity of measurement
        coord: int
            structural frame coordinate
        itype:
            depreciated

        Returns
        -------

        """
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.builders[coord].add_strike_and_dip(pos, strike, dip, polarity)

    def add_tangent_constraint(self, pos, val,coord = None, itype= None):
        # if itype == 'gx':
        #     self.data[0].append(TGPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gy':
        #     self.data[1].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gz':
        #     self.data[2].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        pass

    def add_tangent_constraint_angle(self,pos,s,d,itype):
        # if itype == 'gx':
        #     self.data[0].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gy':
        #     self.data[1].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gz':
        #     self.data[2].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        pass

    def build(self, solver='pyamg', frame=StructuralFrame, **kwargs):
        """
        Build the structural frame
        Parameters
        ----------
        solver solver to use
        frame - type of frame to build StructuralFrame or FoldFrame
        kwargs

        Returns
        -------

        """
        gxxgy = 1
        gxxgz = 1
        gyxgz = 1
        if 'gxxgy' in kwargs:
            gxxgy = kwargs['gxxgy']
        if 'gxxgz' in kwargs:
            gxxgz = kwargs['gxxgz']
        if 'gyxgz' in kwargs:
            gyxgz = kwargs['gyxgz']
        regularisation=kwargs.pop('regularisation',[5.,5.,5.])
        # initialise features as none then where data exists build
        gx_feature = None
        gy_feature = None
        gz_feature = None

        if len(self.builders[0].data) > 0:
            logger.debug("Building structural frame coordinate 0")
            gx_feature = self.builders[0].build(solver=solver,regularisation=regularisation[0],**kwargs)
            # remove fold from kwargs
            fold = kwargs.pop('fold',False)
        if gx_feature is None:
            logger.warning("Not enough constraints for structural frame coordinate 0, \n"
                  "Add some more and try again.")
        if len(self.builders[0].data) > 0:
            logger.debug("Building structural frame coordinate 2")
            # if gy_feature is not None:
            #     self.builders[2].interpolator.add_gradient_orthogonal_constraint(
            #         np.arange(0, self.support.n_elements),
            #         gy_feature.evaluate_gradient(self.support.barycentre),
            #         w=gyxgz)
            if gx_feature is not None:
                self.builders[2].add_orthogonal_feature(gx_feature,gxxgz)
            gz_feature = self.builders[2].build(solver=solver,regularisation=regularisation[1],**kwargs)
        if len(self.builders[0].data) > 0:
            logger.debug("Building structural frame coordinate 1")
            if gx_feature is not None:
                self.builders[1].add_orthogonal_feature(gx_feature,gxxgy)
            if gz_feature is not None:
                self.builders[1].add_orthogonal_feature(gz_feature,gyxgz)
            gy_feature = self.builders[1].build(solver=solver,regularisation=regularisation[2],**kwargs)
        if gy_feature is None:
            logger.warning("Not enough constraints for structural frame coordinate 1, \n"
                  "Add some more and try again.")

        if len(self.builders[2].data) == 0:
            if gy_feature is not None:
                logger.debug("Creating analytical structural frame coordinate 2")
                gz_feature = CrossProductGeologicalFeature(self.name + '_gz',  gy_feature, gx_feature)
            if gy_feature is None or gx_feature is None:
                logger.warning("Not enough constraints for fold frame coordinate 1, \n"
                      "Add some more and try again.")
        # use the frame argument to build a structural frame
        return frame(self.name, [gx_feature, gy_feature, gz_feature])
        
