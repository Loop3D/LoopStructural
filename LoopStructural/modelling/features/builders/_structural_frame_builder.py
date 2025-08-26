"""
structural frame builder
"""

from typing import Union

from LoopStructural.utils.exceptions import LoopException

import numpy as np
import copy

from ....utils import getLogger
from ....datatypes import BoundingBox

logger = getLogger(__name__)


from ....modelling.features.builders import GeologicalFeatureBuilder
from ....modelling.features.builders import FoldedFeatureBuilder
from ....modelling.features import StructuralFrame


class StructuralFrameBuilder:
    def __init__(
        self,
        interpolatortype: Union[str, list],
        bounding_box: BoundingBox,
        nelements: Union[int, list] = 1000,
        frame=StructuralFrame,
        model=None,
        **kwargs,
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
        self.model = model
        # self.region = 'everywhere'
        self.builders = []
        if "name" in kwargs:
            self.name = kwargs["name"]
            kwargs.pop("name")
        self.data = [[], [], []]
        self.fold = kwargs.pop("fold", None)
        # list of interpolators
        # self.interpolators = []
        # Create the interpolation objects by copying the template
        if isinstance(interpolatortype, str):
            interpolatortype = [interpolatortype, interpolatortype, interpolatortype]
        if not isinstance(interpolatortype, list):
            raise LoopException(
                f"interpolatortype is {type(interpolatortype)} and must be either a string or a list of strings"
            )
        if isinstance(nelements, (int, float)):
            nelements = [nelements, nelements, nelements]
        if not isinstance(nelements, list):
            raise LoopException(
                f"nelements is {type(nelements)} and must be either a int or a list of ints"
            )
        # self.builders
        if self.fold:
            self.builders.append(
                FoldedFeatureBuilder(
                    interpolatortype[0],
                    bounding_box,
                    self.fold,
                    nelements=nelements[0],
                    name=f"{self.name}__0",
                    **kwargs,
                )
            )
        else:
            self.builders.append(
                GeologicalFeatureBuilder(
                    interpolatortype[0],
                    bounding_box,
                    nelements[0],
                    name=f"{self.name}__0",
                    **kwargs,
                )
            )  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureBuilder(
                interpolatortype[1],
                bounding_box,
                nelements[1],
                name=f"{self.name}__1",
                **kwargs,
            )
        )  # ,region=self.region))
        self.builders.append(
            GeologicalFeatureBuilder(
                interpolatortype[2],
                bounding_box,
                nelements[2],
                name=f"{self.name}__2",
                **kwargs,
            )
        )  # ,region=self.region))

        self._frame = frame(
            name=self.name,
            features=[
                self.builders[0].feature,
                self.builders[1].feature,
                self.builders[2].feature,
            ],
            fold=self.fold,
            model=self.model,
        )
        self._frame.builder = self

    @classmethod
    def from_feature_builder(cls, feature_builder, **kwargs):
        """
        Create a structural frame builder from an existing feature builder

        Parameters
        ----------
        feature_builder - a geological feature builder
        kwargs

        Returns
        -------

        """
        if not isinstance(feature_builder, GeologicalFeatureBuilder):
            raise LoopException(
                f"feature_builder is {type(feature_builder)} and must be a GeologicalFeatureBuilder"
            )
        if hasattr(feature_builder, 'fold'):
            logger.warning("feature builder has a fold - using this to create a folded frame")
            kwargs['fold'] = copy.deepcopy(feature_builder.fold)
        builder  = cls(
            interpolatortype=[feature_builder.interpolator.type]*3,
            bounding_box=feature_builder.model.bounding_box,
            nelements=[feature_builder.interpolator.n_elements]*3,
            name=feature_builder.name,
            model=feature_builder.model,
            **kwargs
        )
        builder.add_data_from_data_frame(feature_builder.data)
        return builder

    @property
    def build_arguments(self):
        return self.builders[0].build_arguments
    def update_build_arguments(self, kwargs):
        for i in range(3):
            self.builders[i].update_build_arguments(kwargs)
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
            self.builders[i].add_data_from_data_frame(data_frame.loc[data_frame["coord"] == i, :])

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
        if isinstance(regularisation, (int, float)):
            regularisation = np.zeros(3) + regularisation
        logger.info(f"Setting regularisation to {regularisation}")

        # initialise features as none then where data exists build
        if len(self.builders[0].data) > 0:
            logger.info(f"Building {self.name} coordinate 0")
            kwargs["regularisation"] = regularisation[0]
            self.builders[0].update_build_arguments(kwargs)
            kwargs.pop("fold", None)

        # make sure that all of the coordinates are using the same region
        if len(self.builders[2].data) > 0:
            logger.info(f"Building {self.name} coordinate 2")
            if w2 > 0:
                self.builders[2].add_orthogonal_feature(self.builders[0].feature, w2, step=step)
            kwargs["regularisation"] = regularisation[2]
            self.builders[2].update_build_arguments(kwargs)

        if len(self.builders[1].data) > 0:
            logger.info(f"Building {self.name} coordinate 1")
            if w1 > 0:
                self.builders[1].add_orthogonal_feature(self.builders[0].feature, w1, step=step)
            if w3 > 0 and len(self.builders[2].data) > 0:
                self.builders[1].add_orthogonal_feature(self.builders[2].feature, w2, step=step)
            kwargs["regularisation"] = regularisation[1]
            self.builders[1].update_build_arguments(kwargs)

        if len(self.builders[2].data) == 0:
            from LoopStructural.modelling.features import (
                CrossProductGeologicalFeature,
            )

            logger.debug("Creating analytical structural frame coordinate 2")
            c3 = CrossProductGeologicalFeature(self.name + "__2", self._frame[0], self._frame[1])
            self._frame[2] = c3

        # use the frame argument to build a structural frame

    def update(self):
        for i in range(3):
            self.builders[i].update()

    def up_to_date(self, callback=None):
        for i in range(3):
            self.builders[i].up_to_date(callback=callback)

    def set_not_up_to_date(self, caller):

        for i in range(3):
            self.builders[i].set_not_up_to_date(caller)
