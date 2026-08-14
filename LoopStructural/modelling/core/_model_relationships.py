"""Fault/unconformity relationship bookkeeping for GeologicalModel (see API.md).

Extracted from GeologicalModel to separate feature-stack relationship logic
(domain faults, unconformities) from feature-container state and
construction. GeologicalModel's @public_api methods (``add_unconformity``,
``add_onlap_unconformity``) stay defined directly on the class -- their
__qualname__ is part of the CI-checked stable API surface -- and delegate to
the staticmethods here. The private ``_add_faults``/``_add_domain_fault_*``
helpers also delegate, since they're called throughout GeologicalModel's
feature-construction methods.
"""

from ...modelling.features import FeatureType, UnconformityFeature
from ...utils import getLogger

logger = getLogger(__name__)


class FeatureRelationshipManager:
    @staticmethod
    def add_faults(model, feature_builder, features=None):
        """Adds all existing faults to a geological feature builder

        Parameters
        ----------
        model : GeologicalModel
        feature_builder : GeologicalFeatureBuilder/StructuralFrameBuilder
            The feature buider to add the faults to
        features : list, optional
            A specific list of features rather than all features in the model
        """
        if features is None:
            features = model.features
        for f in reversed(features):
            if isinstance(f, str):
                f = model.__getitem__(f)
            if f.type == FeatureType.FAULT:
                feature_builder.add_fault(f)

    @staticmethod
    def add_domain_fault_above(model, feature):
        """
        Looks through the feature list and adds any domain faults to the feature. The domain fault masks everything
        where the fault scalar field is < 0 as being active when added to feature.

        Parameters
        ----------
        model : GeologicalModel
        feature : GeologicalFeatureBuilder
            the feature being added to the model where domain faults should be added
        """
        for f in reversed(model.features):
            if f.name == feature.name:
                continue
            if f.type == "domain_fault":
                feature.add_region(lambda pos, fault=f: fault.evaluate_value(pos) < 0)
                break

    @staticmethod
    def add_domain_fault_below(model, domain_fault):
        """
        Looks through the feature list and adds any the domain_fault to the features
        that already exist in the stack until an unconformity is reached. domain faults
        to the feature. The domain fault masks everything where the fault scalar field
         is < 0 as being active when added to feature.

        Parameters
        ----------
        model : GeologicalModel
        domain_fault : GeologicalFeatureBuilder
            the feature being added to the model where domain faults should be added
        """
        for f in reversed(model.features):
            if f.name == domain_fault.name:
                continue
            f.add_region(lambda pos: domain_fault.evaluate_value(pos) > 0)
            if f.type == FeatureType.UNCONFORMITY:
                break

    @staticmethod
    def add_unconformity_above(model, feature):
        """
        Adds a region to the feature to prevent the value from being
        interpolated where the unconformities exists above e.g.
        if there is another feature above and the unconformity is at 0
        then the features added below (after) will only be visible where the
        uncomformity is <0

        Parameters
        ----------
        model : GeologicalModel
        feature - GeologicalFeature
        """

        if feature.type == FeatureType.FAULT:
            return
        for f in reversed(model.features):
            if f.type == FeatureType.UNCONFORMITY and f.name != feature.name:
                logger.info(f"Adding {f.name} as unconformity to {feature.name}")
                feature.add_region(f)
            if f.type == FeatureType.ONLAPUNCONFORMITY and f.name != feature.name:
                feature.add_region(f)
                break

    @staticmethod
    def add_unconformity(model, feature, value, index=None):
        """
        Use an existing feature to add an unconformity to the model.

        Parameters
        ----------
        model : GeologicalModel
        feature : GeologicalFeature
            existing geological feature
        value : float
            scalar value of isosurface that represents

        Returns
        -------
        unconformity : GeologicalFeature
            unconformity feature
        """
        logger.debug(f"Adding {feature.name} as unconformity at {value}")
        if feature is None:
            logger.warning("Cannot add unconformtiy, base feature is None")
            return
        # look backwards through features and add the unconformity as a region until
        # we get to an unconformity
        uc_feature = UnconformityFeature(feature, value)
        feature.add_region(uc_feature.inverse())
        for f in reversed(model.features):
            if f.type == FeatureType.UNCONFORMITY:
                logger.debug(f"Reached unconformity {f.name}")
                break
            logger.debug(f"Adding {uc_feature.name} as unconformity to {f.name}")
            if f.type == FeatureType.FAULT or f.type == FeatureType.INACTIVEFAULT:
                continue
            if f == feature:
                continue
            else:
                f.add_region(uc_feature)
        # now add the unconformity to the feature list
        model._add_feature(uc_feature, index=index)
        return uc_feature

    @staticmethod
    def add_onlap_unconformity(model, feature, value, index=None):
        """
        Use an existing feature to add an unconformity to the model.

        Parameters
        ----------
        model : GeologicalModel
        feature : GeologicalFeature
            existing geological feature
        value : float
            scalar value of isosurface that represents

        Returns
        -------
        unconformity_feature : GeologicalFeature
            the created unconformity
        """
        feature.regions = []
        uc_feature = UnconformityFeature(feature, value, False, onlap=True)
        feature.add_region(uc_feature.inverse())
        for f in reversed(model.features):
            if f.type in (FeatureType.UNCONFORMITY, FeatureType.ONLAPUNCONFORMITY):
                logger.debug(f"Reached unconformity {f.name}")
                break
            if f.type == FeatureType.FAULT or f.type == FeatureType.INACTIVEFAULT:
                continue
            if f != feature:
                f.add_region(uc_feature)
        model._add_feature(uc_feature.inverse(), index=index)

        return uc_feature
