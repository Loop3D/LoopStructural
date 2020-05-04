import logging

import numpy as np

logger = logging.getLogger(__name__)


class GeologicalFeature:
    def __init__(self, name, interpolator, builder=None, data=None, region=None, type=None, faults=[]):
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
        self.interpolator = interpolator
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
        self.model = None
    def __str__(self):
        return self.name

    def set_model(self, model):
        self.model = model

    def toggle_faults(self):
        """
        Turn the fault off for a feature
        This function is only really used for debugging or creating methods
        explanation figures

        Returns
        -------

        """
        self.faults_enabled = ~self.faults_enabled

    def add_region(self, region):
        """
        Adds a region where the geological feature is active to the model.

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
        Evaluate the scalar field value of the geological feature at the locations
        specified

        Parameters
        ----------
        evaluation_points : numpy array
            location to evaluate the scalar value

        Returns
        -------
        values : numpy array
            numpy array containing evaluated values

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
        v[mask] = self.interpolator.evaluate_value(evaluation_points[mask, :])
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
        v[mask, :] = self.interpolator.evaluate_gradient(evaluation_points)

        return v

    def evaluate_gradient_misfit(self):
        grad = self.interpolator.get_gradient_constraints()
        norm = self.interpolator.get_norm_constraints()
        dot = []
        print(grad.shape)
        print(norm.shape)
        if grad.shape[0] > 0:
            model_grad = self.evaluate_gradient(grad[:,:3])
            dot.append(np.einsum('ij,ij->i',model_grad,grad[:,:3:6]).tolist())

        if norm.shape[0] > 0:
            print(norm[:,:3])
            model_norm = self.evaluate_gradient(norm[:, :3])
            print(model_norm)
            dot.append(np.einsum('ij,ij->i', model_norm, norm[:,:3:6]))

        return np.array(dot)

    def evaluate_value_misfit(self):
        locations = self.interpolator.get_value_constraints()
        return locations[:,3] - self.evaluate_value(locations[:,:3])

    def mean(self):
        """
        Calculate average of the support values

        Returns
        -------

        """
        if self.model is None:
            return 0
        return np.mean(self.evaluate_value(self.model.regular_grid((10,10,10))))
        # return np.nanmean(self.scalar_field.get_node_values())

    def min(self):
        """

        Returns
        -------

        """
        if self.model is None:
            return 0
        return np.min(
            self.evaluate_value(self.model.regular_grid((10, 10, 10))))
        #       return np.nanmin(self.scalar_field.get_node_values())

    def max(self):
        """
        Calculate average of the support values

        Returns
        -------

        """
        if self.model is None:
            return 0
        return np.max(
            self.evaluate_value(self.model.regular_grid((10, 10, 10))))
        #return np.nanmax(self.scalar_field.get_node_values())

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
        self.interpolator.up_to_date = False
        self.interpolator.update()

    def get_interpolator(self):
        """
        Get the interpolator used to build this feature

        Returns
        -------
        GeologicalInterpolator
        """
        return self.interpolator

    # def get_node_values(self):
    #     """
    #     Get the node values of the support used to build this interpolator
    #     if the
    #     interpolator is a discrete interpolator
    #
    #     Returns
    #     -------
    #     numpy array of values
    #     """
    #     # if interpolator.type =
    #     return self.scalar_field.get_node_values()
