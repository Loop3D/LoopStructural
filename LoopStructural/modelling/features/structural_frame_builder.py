"""
structural frame builder
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

from LoopStructural.modelling.features.cross_product_geological_feature \
    import \
    CrossProductGeologicalFeature
from LoopStructural.modelling.features import \
    GeologicalFeatureInterpolator
from LoopStructural.modelling.features import StructuralFrame


class StructuralFrameBuilder:
    """[summary]

    [extended_summary]
    """
    def __init__(self, interpolator=None, interpolators=None, **kwargs):
        """
        Class for building a structural frame - has functions to set up the
        interpolator with
        data and also orthogonality constraints. Can build a generic
        structural frame or a
        subclass of the structural frame if the kwarg `frame` is specified

        Parameters
        ----------
        interpolator - a template interpolator for the frame
        kwargs
        """
        self.support = None
        self.fault_event = None
        self.name = 'Undefined'
        # self.region = 'everywhere'
        self.builders = []
        if 'name' in kwargs:
            self.name = kwargs['name']
            kwargs.pop('name')
        self.data = [[], [], []]
        # list of interpolators
        # self.interpolators = []
        # Create the interpolation objects by copying the template
        if interpolators is None:
            if interpolator is not None:
                interpolators = []
                interpolators.append(interpolator)
                interpolators.append(interpolator.copy())
                interpolators.append(interpolator.copy())
            else:
                raise BaseException("Missing interpolator")
        # self.builders
        self.builders.append(
            GeologicalFeatureInterpolator(interpolators[0],
                                          name=self.name + '_0',
                                          **kwargs))  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureInterpolator(interpolators[1],
                                          name=self.name + '_1',
                                          **kwargs))  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureInterpolator(interpolators[2],
                                          name=self.name + '_2',
                                          **kwargs))  # ,region=self.region))

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
        for i in range(3):
            self.builders[i].add_data_from_data_frame(data_frame.loc[data_frame['coord'] == i,:])


    def build(self, w1=1., w2=1., w3=1., frame=StructuralFrame, **kwargs):
        """
        Build the structural frame
        Parameters
        ----------
        solver solver to use
        frame - type of frame to build StructuralFrame or FoldFrame
        w3
        w2
        w1
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
        # set regularisation so the the main surface (foliation, fault) is smooth
        # and the fields are allowed to vary more
        regularisation = kwargs.pop('regularisation', [1., 1., 1.])
        # initialise features as none then where data exists build
        gx_feature = None
        gy_feature = None
        gz_feature = None
        fold = None
        if len(self.builders[0].data) > 0:
            logger.info("Building %s coordinate 0"%self.name)
            gx_feature = self.builders[0].build(regularisation=regularisation[
                                                    0], **kwargs)
            # remove fold from kwargs

            fold = kwargs.pop('fold', False)
        if gx_feature is None:
            logger.warning(
                "Not enough constraints for structural frame coordinate 0, \n"
                "Add some more and try again.")
        # make sure that all of the coordinates are using the same region
        if gx_feature is not None and gx_feature.get_interpolator().region_function is not None:
            self.builders[1].interpolator.set_region(gx_feature.get_interpolator().region_function)
            self.builders[2].interpolator.set_region(gx_feature.get_interpolator().region_function)
            if 'data_region' in kwargs:
                kwargs.pop('data_region')
            if 'region' in kwargs:
                kwargs.pop('region')
        if len(self.builders[2].data) > 0:
            logger.info("Building %s coordinate 2"%self.name)
            # if gy_feature is not None:
            #     self.builders[
            #     2].interpolator.add_gradient_orthogonal_constraint(
            #         np.arange(0, self.support.n_elements),
            #         gy_feature.evaluate_gradient(self.support.barycentre),
            #         w=gyxgz)
            if gx_feature is not None:
                self.builders[2].add_orthogonal_feature(gx_feature, gxxgz)
            gz_feature = self.builders[2].build(regularisation=regularisation[2], **kwargs)

        if len(self.builders[1].data) > 0:
            logger.info("Building %s coordinate 1"%self.name)
            if gx_feature is not None:
                self.builders[1].add_orthogonal_feature(gx_feature, gxxgy)
            if gz_feature is not None:
                self.builders[1].add_orthogonal_feature(gz_feature, gyxgz)
            gy_feature = self.builders[1].build(regularisation=regularisation[1], **kwargs)

        if gy_feature is None:
            logger.warning(
                "Not enough constraints for structural frame coordinate 1, \n"
                "Add some more and try again.")

        if len(self.builders[2].data) == 0:
            if gy_feature is not None:
                logger.debug(
                    "Creating analytical structural frame coordinate 2")
                gz_feature = CrossProductGeologicalFeature(self.name + '_2',
                                                           gy_feature,
                                                           gx_feature)
            if gy_feature is None or gx_feature is None:
                logger.warning(
                    "Not enough constraints for fold frame coordinate 1, \n"
                    "Add some more and try again.")
        # use the frame argument to build a structural frame
        return frame(self.name, [gx_feature, gy_feature, gz_feature],fold=fold)
