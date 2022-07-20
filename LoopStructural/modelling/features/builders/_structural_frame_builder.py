"""
structural frame builder
"""

import logging

import numpy as np

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


from LoopStructural.modelling.features.builders import GeologicalFeatureBuilder
from LoopStructural.modelling.features.builders import FoldedFeatureBuilder
from LoopStructural.modelling.features import StructuralFrame


class StructuralFrameBuilder:
    def __init__(
        self, interpolator=None, interpolators=None, frame=StructuralFrame, **kwargs
    ):
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
        self.name = "Undefined"
        # self.region = 'everywhere'
        self.builders = []
        if "name" in kwargs:
            self.name = kwargs["name"]
            kwargs.pop("name")
        self.data = [[], [], []]
        self.fold = kwargs.get("fold", None)
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
        if "fold" in kwargs:
            self.builders.append(
                FoldedFeatureBuilder(interpolators[0], name=f"{self.name}__0", **kwargs)
            )
        else:
            self.builders.append(
                GeologicalFeatureBuilder(
                    interpolators[0], name=f"{self.name}__0", **kwargs
                )
            )  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureBuilder(interpolators[1], name=f"{self.name}__1", **kwargs)
        )  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureBuilder(interpolators[2], name=f"{self.name}__2", **kwargs)
        )  # ,region=self.region))

        self._frame = frame(
            self.name,
            [
                self.builders[0].feature,
                self.builders[1].feature,
                self.builders[2].feature,
            ],
            fold=kwargs.get("fold", None),
        )
        self._frame.builder = self

    @property
    def frame(self):
        return self._frame

    def __getitem__(self, item):
        return self.builders[item]

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
        for i in range(3):
            self.builders[i].add_fault(fault)

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
            self.builders[i].add_data_from_data_frame(
                data_frame.loc[data_frame["coord"] == i, :]
            )

    def setup(self, w1=1.0, w2=1.0, w3=1.0, **kwargs):
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
        step = kwargs.get("step", 10)
        if "gxxgy" in kwargs:
            logger.warning("gxxgy deprecated please use w1")
            w1 = kwargs["gxxgy"]
        if "gxxgz" in kwargs:
            logger.warning("gxxgz deprecated please use w2")
            w2 = kwargs["gxxgz"]
        if "gyxgz" in kwargs:
            logger.warning("gyxgz deprecated please use w3")
            w3 = kwargs["gyxgz"]

        # set regularisation so the the main surface (foliation, fault) is smooth
        # and the fields are allowed to vary more
        regularisation = kwargs.pop("regularisation", [1.0, 1.0, 1.0])
        # initialise features as none then where data exists build
        if len(self.builders[0].data) > 0:
            logger.info(f"Building {self.name} coordinate 0")
            kwargs["regularisation"] = regularisation[0]
            self.builders[0].build_arguments = kwargs
            fold = kwargs.pop("fold", None)

        # make sure that all of the coordinates are using the same region
        if len(self.builders[2].data) > 0:
            logger.info(f"Building {self.name} coordinate 2")
            if w2 > 0:
                self.builders[2].add_orthogonal_feature(
                    self.builders[0].feature, w2, step=step
                )
            kwargs["regularisation"] = regularisation[2]
            self.builders[2].build_arguments = kwargs

        if len(self.builders[1].data) > 0:
            logger.info(f"Building {self.name} coordinate 1")
            if w1 > 0:
                self.builders[1].add_orthogonal_feature(
                    self.builders[0].feature, w1, step=step
                )
            if w3 > 0 and len(self.builders[2].data) > 0:
                self.builders[1].add_orthogonal_feature(
                    self.builders[2].feature, w2, step=step
                )
            kwargs["regularisation"] = regularisation[1]
            self.builders[1].build_arguments = kwargs

        if len(self.builders[2].data) == 0:
            from LoopStructural.modelling.features import (
                CrossProductGeologicalFeature,
            )

            logger.debug("Creating analytical structural frame coordinate 2")
            c3 = CrossProductGeologicalFeature(
                self.name + "_2", self._frame[0], self._frame[1]
            )
            self._frame[2] = c3

        # use the frame argument to build a structural frame

    def update(self):
        for i in range(3):
            self.builders[i].update()

    def up_to_date(self, callback=None):
        for i in range(3):
            self.builders[i].up_to_date(callback=callback)
