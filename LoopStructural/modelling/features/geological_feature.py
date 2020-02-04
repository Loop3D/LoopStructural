import logging

import numpy as np

logger = logging.getLogger(__name__)


class GeologicalFeature:
    def __init__(self, name, support, builder=None, data=None, region=None, type=None, faults=[]):
        """
        Geological feature is class that is used to represent a geometrical element in a geological
        model. For example foliations, fault planes, fold rotation angles etc. The feature has a support
        which

        Parameters
        ----------
        name: string
        support
        builder
        data
        region

        Attributes
        ----------
        name : string
            should be a unique name for the geological feature
        support : a ScalarField
            holds the property values for the feature and links to the
            support geometry
        data : list
            list containing geological data
        region : list of boolean functions defining whether the feature is
            active
        faults : list of faults that affect this feature
        """
        self.name = name
        self.support = support
        self.ndim = 1
        self.data = data
        self.builder = builder
        self.region = region
        self.regions = []
        self.type = type
        self.faults = faults
        self.faults_enabled = True
        if region is None:
            self.region = 'everywhere'

    def __str__(self):
        return self.name

    def toggle_faults(self):
        self.faults_enabled = ~self.faults_enabled

    def add_region(self, region):
        """

        Parameters
        ----------
        region : boolean function(x,y,z)
                - returns true if inside region, false if outside
                can be passed as a lambda function e.g.
                lambda pos : feature.evaluate_value(pos) > 0

        Returns
        -------

        """
        self.regions.append(region)

    def set_builder(self, builder):
        """

        Parameters
        ----------
        builder

        Returns
        -------

        """
        self.builder = builder

    def evaluate_value(self, evaluation_points):
        """

        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """

        # check if the points are within the display region
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask, r(evaluation_points))
        # apply faulting after working out which regions are visible
        if self.faults_enabled:
            for f in self.faults:
                evaluation_points = f.apply_to_points(evaluation_points)
        v[mask] = self.support.evaluate_value(evaluation_points[mask, :])
        return v

    def evaluate_gradient(self, evaluation_points):
        """

        Parameters
        ----------
        locations : numpy array
            location where the gradient is being evaluated

        Returns
        -------

        """
        v = np.zeros(evaluation_points.shape)
        v[:] = np.nan
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask, r(evaluation_points))

        # apply faulting after working out which regions are visible
        if self.faults_enabled:
            for f in self.faults:
                evaluation_points = f.apply_to_points(evaluation_points)
        v[mask, :] = self.support.evaluate_gradient(evaluation_points)

        return v

    def mean(self):
        """
        Calculate average of the support values

        Returns
        -------

        """
        return np.nanmean(self.support.get_node_values())

    def min(self):
        """

        Returns
        -------

        """
        return np.nanmin(self.support.get_node_values())

    def max(self):
        """
        Calculate average of the support values

        Returns
        -------

        """
        return np.nanmax(self.support.get_node_values())

    def update(self):
        """
        Calculate average of the support values

        Returns
        -------

        """
        # re-run the interpolator and update the support.
        # this is a bit clumsy and not abstract, i think
        # if evaluating the property doesn't require the dictionary on
        # the nodes and actually just uses the interpolator values this
        # would be
        # much better.
        self.support.interpolator.up_to_date = False
        self.support.interpolator.update()
        self.support.update_property(self.support.interpolator.c)

    def get_interpolator(self):
        """
        Get the interpolator used to build this feature

        Returns
        -------
        GeologicalInterpolator
        """
        return self.support.interpolator

    def get_node_values(self):
        """
        Get the node values of the support used to build this interpolator
        if the
        interpolator is a discrete interpolator

        Returns
        -------
        numpy array of values
        """
        return self.support.get_node_values()

    def slice(self, **kwargs):
        logger.error("function has been removed, please use the modelviewer class")
        return
