"""
Feature builder
"""
import copy
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from LoopStructural.utils.helper import xyz_names, val_name, normal_vec_names, \
    weight_name, gradient_vec_names, tangent_vec_names, interface_name
from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.utils.helper import get_data_bounding_box_map as get_data_bounding_box


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
        self.interpolator = interpolator
        self.name = name
        self.interpolator.set_property_name(self.name)
        # everywhere region is just a lambda that returns true for all locations
        if region is None:
            self.region = lambda pos: np.ones(pos.shape[0], dtype=bool)
        else:
            self.region = region
        header = xyz_names()+val_name()+gradient_vec_names()+\
                 normal_vec_names()+tangent_vec_names()+weight_name()
        self.data = pd.DataFrame(columns=header)
        self.faults = []
        self.data_added = False
        self.interpolator.set_region(region=self.region)

    def update(self):
        pass

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
        self.faults.append(fault)

    def add_data_from_data_frame(self, data_frame):
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

    def add_orthogonal_feature(self, feature, w=1., region=None):
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

        Returns
        -------

        """
        self.interpolator.add_gradient_orthogonal_constraint(
            self.interpolator.support.barycentre(),
            feature.evaluate_gradient(self.interpolator.support.barycentre()),
            w=w
        )

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
            return np.zeros(0, 7)

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
            return np.zeros(0, 7)

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
            return np.zeros(0,7)

    def get_interface_constraints(self):
        mask = np.all(~np.isnan(self.data.loc[:, interface_name()].to_numpy()), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[
                mask, xyz_names() + interface_name() + weight_name()].to_numpy()
        else:
            return np.zeros(0,5)

    def get_data_locations(self):
        """
        Get only the location for all data points

        Returns
        -------

        """
        return self.data.loc[:, xyz_names()].to_numpy()

    def build(self, fold=None, fold_weights=None, data_region=None, **kwargs):
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

        if data_region is not None:
            xyz = self.get_data_locations()
            bb, region = get_data_bounding_box(xyz, data_region)
            self.interpolator.set_region(region=region)
        if not self.data_added:
            self.add_data_to_interpolator(**kwargs)

        # moving this to init because it needs to be done before constraints
        # are added?
        if fold is not None:
            logger.info("Adding fold to %s" % self.name)
            self.interpolator.fold = fold
            # if we have fold weights use those, otherwise just use default
            if fold_weights is None:
                self.interpolator.add_fold_constraints()
            else:
                self.interpolator.add_fold_constraints(fold_weights)
            if 'cgw' not in kwargs:
                kwargs['cgw'] = 0.

        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(**kwargs)
        return GeologicalFeature(self.name,
                                 self.interpolator,
                                 builder=self, data=self.data,
                                 region=self.region,
                                 faults=self.faults,
                                 fold = fold
                                 )
