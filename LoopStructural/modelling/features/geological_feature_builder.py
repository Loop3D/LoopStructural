"""
Feature builder
"""
import copy
import logging

import numpy as np
import pandas as pd

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

from LoopStructural.utils.helper import xyz_names, val_name, normal_vec_names, \
    weight_name, gradient_vec_names, tangent_vec_names, interface_name
from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.utils.helper import get_data_bounding_box_map as get_data_bounding_box
from LoopStructural.utils import get_data_axis_aligned_bounding_box
from LoopStructural.utils import RegionEverywhere

class GeologicalFeatureInterpolator:
    """[summary]

    [extended_summary]
    """
    def __init__(self, interpolator, name='Feature', region=None, **kwargs):
        """
        Constructor for a GeologicalFeatureInterpolator

        Parameters
        ----------
        interpolator : GeologicalInterpolator
            An empty GeologicalInterpolator
        region : lambda function
            defining whether the location (xyz) should be included in the
        kwargs - name of the feature, region to interpolate the feature
        """
        self._interpolator = interpolator
        self._name = name
        self._interpolator.set_property_name(self._name)
        # everywhere region is just a lambda that returns true for all locations
        if region is None:
            self.region = RegionEverywhere()
        else:
            self.region = region
        header = xyz_names()+val_name()+gradient_vec_names()+\
                 normal_vec_names()+tangent_vec_names()+weight_name()
        self.data = pd.DataFrame(columns=header)
        self.faults = []
        self.data_added = False
        self._interpolator.set_region(region=self.region)
        self._feature = None
        self._up_to_date = False
        self._build_arguments = {}
        self.fold = None
        self._feature = GeologicalFeature(self._name,
                                 self._interpolator,
                                 builder=self, 
                                 region=self.region,
                                 faults=self.faults,
                                 fold = self.fold
                                 )
        self._orthogonal_features = {}
        self._equality_constraints = {}
    @property
    def feature(self):
        return self._feature

    @property
    def build_arguments(self):
        return self._build_arguments 
    
    @build_arguments.setter
    def build_arguments(self, build_arguments):
        # self._build_arguments = {}
        for k, i in build_arguments.items():
            self._build_arguments[k] = i

    def update(self):
        self.build(**self.build_arguments)
    @property
    def name(self):
        return self._name
    @property
    def interpolator(self):
        return self._interpolator

    def up_to_date(self):
        #has anything changed in the builder since we built the feature? if so update
        if self._up_to_date == False:
            self.update()
        #check if the interpolator is up to date, if not solve
        if self._interpolator.up_to_date == False:
            self.update()
        
        
    def add_fault(self, fault):
        """
        Add a fault to the geological feature builder

        Parameters
        ----------
        fault : FaultSegment
            A faultsegment to add to the geological feature

        Returns
        -------

        """
        self._up_to_date = False
        self.faults.append(fault)

    def add_data_from_data_frame(self, data_frame, overwrite = False):
        """
        Extract data from a pandas dataframe with columns for

        Parameters
        ----------
        data_frame : pd.DataFrame
            a dataframe containing the data to be added

        Returns
        -------

        """
        self.data = data_frame.copy()

    def add_orthogonal_feature(self, feature, w=1., region=None,step=1,B=0):
        """
        Add a constraint to the interpolator so that the gradient of an exisitng feature is orthogonal
        to the feature being built. E.g. dot product between gradients should be = 0

        Parameters
        ----------
        feature : GeologicalFeature
            feature which we want to be orthogonal to
        w :  double
            how much to weight in least squares sense
        region : unused
        step : int
            numpy slicing step size to see how many tetras to add
        
        Returns
        -------

        Notes
        -----
        The constraint can be applied to a random subset of the tetrahedral elements in the mesh
        in theory this shu
        """
        try:
            step = int(step) #cast as int in case it was a float
        except ValueError:
            logger.error("Cannot cast {} as integer, setting step to 1".format(step))
            step = 1
        self._orthogonal_features[feature.name] = [feature,w,region,step,B]
        self._up_to_date = False

        
    
    def add_data_to_interpolator(self, constrained=False, force_constrained=False, **kwargs):
        """
        Iterates through the list of data and applies any faults active on the
        data in the order they are added

        Parameters
        -----------
        constrained : boolean
        force_constrained : boolean

        Returns
        -------

        """
        # first move the data for the fault
        logger.info("Adding %i faults to %s" % (len(self.faults), self.name))
        data = self.data.copy()
        # convert data locations to numpy array and then update
        for f in self.faults:
            data.loc[:,xyz_names()] = f.apply_to_points(
                self.get_data_locations())
        # Now check whether there are enough constraints for the
        # interpolator to be able to solve
        # we need at least 2 different value points or a single norm
        # constraint. If there are not enough
        # try converting grad to norms, if still not enough send user an error
        if constrained:
            # Change normals to gradients
            mask = np.all(~np.isnan(data.loc[:, normal_vec_names()]),axis=1)
            if mask.shape[0] > 0:
                data.loc[mask, gradient_vec_names()] = data.loc[mask,
                                                            normal_vec_names()].to_numpy()
                data.loc[mask, normal_vec_names()] = np.nan
        if self.get_norm_constraints().shape[0] > 0:
            constrained = True

        if np.unique(self.get_value_constraints()[:,3]).shape[0]>1:
            constrained = True

        if not constrained or force_constrained:
            # change gradient constraints to normal vector constraints
            mask = np.all(~np.isnan(data.loc[:, gradient_vec_names()]), axis=1)
            if mask.shape[0] > 0:

                data.loc[mask, normal_vec_names()] = data.loc[mask,
                                                            gradient_vec_names()].to_numpy()
                data.loc[mask, gradient_vec_names()] = np.nan
                logger.info(
                    "Setting gradient points to norm constraints")
                constrained = True
                mask = np.all(
                    ~np.isnan(data.loc[:, normal_vec_names()].to_numpy()),
                    axis=1)

        if not constrained:
            logger.error("Not enough constraints for scalar field add more")
        # self.interpolator.reset()
        mask = ~np.isnan(data.loc[:,val_name()].to_numpy())
        # add value constraints
        if mask.shape[0]>0:
            value_data = data.loc[mask[:,0],xyz_names()+val_name()+weight_name()].to_numpy()
            self.interpolator.set_value_constraints(value_data)

        # add gradient constraints
        mask = np.all(~np.isnan(data.loc[:, gradient_vec_names()].to_numpy()), axis=1)
        if mask.shape[0]>0:
            gradient_data = data.loc[
            mask, xyz_names() + gradient_vec_names() + weight_name()].to_numpy()
            self.interpolator.set_gradient_constraints(gradient_data)

        # add normal vector data
        mask = np.all(~np.isnan(data.loc[:, normal_vec_names()].to_numpy()), axis=1)
        if mask.shape[0]>0:
            normal_data = data.loc[
                mask, xyz_names() + normal_vec_names() + weight_name()].to_numpy()
            self.interpolator.set_normal_constraints(normal_data)

        # add tangent data
        mask = np.all(~np.isnan(data.loc[:, tangent_vec_names()].to_numpy()), axis=1)
        if mask.shape[0]>0:
            tangent_data = data.loc[
                mask, xyz_names() + tangent_vec_names() + weight_name()].to_numpy()
            self.interpolator.set_tangent_constraints(tangent_data)

        # add interface constraints
        mask = np.all(~np.isnan(data.loc[:, interface_name()].to_numpy()), axis=1)
        if mask.shape[0] > 0:
            interface_data = data.loc[
                mask, xyz_names() + interface_name() + weight_name()].to_numpy()
            self.interpolator.set_interface_constraints(interface_data)

        self.data_added = True
        self._up_to_date = False

    
    def install_gradient_constraint(self):
        for g in self._orthogonal_features.values():
            feature,w,region,step,B = g
            vector = feature.evaluate_gradient(self.interpolator.support.barycentre())
            vector /= np.linalg.norm(vector,axis=1)[:,None]
            element_idx = np.arange(self.interpolator.support.n_elements)
            np.random.shuffle(element_idx)
            self.interpolator.add_gradient_orthogonal_constraint(
                self.interpolator.support.barycentre()[element_idx[::step],:],
                vector[element_idx[::step],:],
                w=w,
                B=B
            )
    def add_equality_constraints(self,feature, region):

        self._equality_constraints[feature.name] = [feature,region]
        self._up_to_date = False


    def install_equality_constraints(self):
        for e in self._equality_constraints.values():
            try:
                # assume all parts of structural frame have the same support
                support = self.interpolator.support
                
                # work out the values of the nodes where we want hard
                # constraints
                idc = np.arange(0, support.n_nodes)[
                    e[1](support.nodes)]
                val = e[0].evaluate_value(
                    support.nodes[
                    e[1](support.nodes), :])
                mask = ~np.isnan(val)
                self.interpolator.add_equality_constraints(
                    idc[mask], val[mask])
            except:
                logger.error("Could not add equality")


    def get_value_constraints(self):
        """
        Get the value constraints for this geological feature

        Returns
        -------
        np.array((N,4),dtype=double)
        """
        header = xyz_names()+val_name()+weight_name()
        mask = ~np.isnan(self.data.loc[:,val_name()].to_numpy())
        return self.data.loc[mask[:,0],header].to_numpy()

    def get_gradient_constraints(self):
        """
        Get the gradient direction constraints

        Returns
        -------
        numpy array
        """
        mask = np.all(
            ~np.isnan(self.data.loc[:, gradient_vec_names()].to_numpy()),
            axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[
                mask, xyz_names() + gradient_vec_names() + weight_name(
                )].to_numpy()
        else:
            return np.zeros((0, 7))

    def get_tangent_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        mask = np.all(
            ~np.isnan(self.data.loc[:, tangent_vec_names()].to_numpy()),
            axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[
                mask, xyz_names() + tangent_vec_names() + weight_name(
                )].to_numpy()
        else:
            return np.zeros((0, 7))

    def get_norm_constraints(self):
        """
        Get the gradient norm constraints

        Returns
        -------
        numpy array
        """
        mask = np.all(~np.isnan(self.data.loc[:, normal_vec_names()].to_numpy()),
                      axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[
                mask, xyz_names() + normal_vec_names() + weight_name(

                )].to_numpy()
        else:
            return np.zeros((0,7))

    def get_interface_constraints(self):
        mask = np.all(~np.isnan(self.data.loc[:, interface_name()].to_numpy()), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[
                mask, xyz_names() + interface_name() + weight_name()].to_numpy()
        else:
            return np.zeros((0,5))

    def get_data_locations(self):
        """
        Get only the location for all data points

        Returns
        -------

        """
        return self.data.loc[:, xyz_names()].to_numpy()
        
    def set_interpolation_geometry(self,origin,maximum, rotation=None):
        """Update the interpolation support geometry to new bounding box

        Parameters
        ----------
        origin : np.array(3)
            origin vector
        maximum : np.array(3)
            maximum vector
        """
        logger.info("Setting mesh origin: {} {} {} ".format(origin[0],origin[1],origin[2]))
        logger.info("Setting mesh maximum: {} {} {}".format(maximum[0],maximum[1],maximum[2]))
        self.interpolator.support.origin = origin
        self.interpolator.support.maximum = maximum
        self.interpolator.support.rotation_xy = rotation
        self._up_to_date = False

        while self.interpolator.nx < 100:
            self.interpolator.support.step_vector=self.interpolator.support.step_vector*0.9
            
    def build(self, fold=None, fold_weights={}, data_region=None, **kwargs):
        """
        Runs the interpolation and builds the geological feature

        Parameters
        ----------
        fold : FoldEvent
        fold_weights : dict
        data_region : double <1
            If not none adds a region around the data points to the interpolation
            with data_region as a buffer
        kwargs

        Returns
        -------

        """


        self.add_data_to_interpolator(**kwargs)
        if data_region is not None:
            xyz = self.interpolator.get_data_locations()
            bb, region = get_data_bounding_box(xyz, data_region)    
            # if self.model.reuse_supports == False:
            if np.any(np.min(bb[0,:])< self.interpolator.support.origin):
                neworigin = np.min([self.interpolator.support.origin,bb[0,:]],axis=0)
                logger.info("Changing origin of support for {} from {} {} {} to {} {} {}"\
                                .format(self.name,self.interpolator.support.origin[0],\
                                    self.interpolator.support.origin[1],self.interpolator.support.origin[2],\
                                                neworigin[0],neworigin[1],neworigin[2]))
                self.interpolator.support.origin = neworigin
            if np.any(np.max(bb[0,:])< self.interpolator.support.maximum):
                newmax = np.max([self.interpolator.support.maximum,bb[0,:]],axis=0)
                logger.info("Changing origin of support for {} from {} {} {} to {} {} {}"\
                        .format(self.name,self.interpolator.support.maximum[0],
                        self.interpolator.support.maximum[1],self.interpolator.support.maximum[2],\
                                                newmax[0],newmax[1],newmax[2]))

                self.interpolator.support.maximum = newmax
            self.interpolator.set_region(region=region)
        # moving this to init because it needs to be done before constraints
        # are added?
        if fold is not None:
            logger.info("Adding fold to %s" % self.name)
            self.interpolator.fold = fold
            # if we have fold weights use those, otherwise just use default
            self.interpolator.add_fold_constraints(**fold_weights)
            if 'cgw' not in kwargs:
                # try adding very small cg
                kwargs['cgw'] = 0.0
        self.install_gradient_constraint()
        self.install_equality_constraints()
        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(**kwargs)
        self._up_to_date = True
        return self._feature
