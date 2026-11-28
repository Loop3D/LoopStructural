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

from __future__ import annotations

from typing import TYPE_CHECKING

from ...modelling.features import FeatureType, GeologicalFeature, UnconformityFeature
from ...modelling.features.builders import GeologicalFeatureBuilder, StructuralFrameBuilder

if TYPE_CHECKING:
    from .geological_model import GeologicalModel
from ...utils import getLogger

logger = getLogger(__name__)


class FeatureRelationshipManager:
    @staticmethod
    def add_faults(
        model: GeologicalModel,
        feature_builder: GeologicalFeatureBuilder | StructuralFrameBuilder,
        features: list[GeologicalFeature | str] | None = None,
    ) -> None:
        """Add all existing faults to a geological feature builder.

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
    def add_domain_fault_above(model: GeologicalModel, feature: GeologicalFeatureBuilder) -> None:
        """Add the active domain fault above the feature.

        The domain fault masks everything
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
            if f.type == FeatureType.DOMAINFAULT:
                # A proper Feature-like region (name/parent), not a bare
                # lambda -- code that walks `self.regions` expecting
                # Feature-like objects (e.g. BaseFeature.surfaces) would
                # otherwise raise AttributeError on this region and
                # silently produce no surfaces for anything a domain fault
                # crops.
                feature.add_region(UnconformityFeature(f, 0, sign=True))
                break

    @staticmethod
    def add_domain_fault_below(model: GeologicalModel, domain_fault: GeologicalFeature) -> None:
        """Add the domain fault to features below it.

        The domain fault is added to features
        that already exist in the stack until an unconformity or another domain
        fault is reached (that boundary is not cropped). domain faults
        to the feature. The domain fault masks everything where the fault scalar field
         is < 0 as being active when added to feature.

        Parameters
        ----------
        model : GeologicalModel
        domain_fault : GeologicalFeatureBuilder
            the feature being added to the model where domain faults should be added
        """
        region = UnconformityFeature(domain_fault, 0, sign=False)
        for f in reversed(model.features):
            if f.name == domain_fault.name:
                continue
            if f.type == FeatureType.UNCONFORMITY or f.type == FeatureType.DOMAINFAULT:
                # An unconformity or domain fault added before this one is
                # the boundary above the features this domain fault crops.
                # It is a stopping point and is not itself cropped, else the
                # boundary surface is removed on one side of this fault.
                # This matches add_unconformity, which also stops before it
                # crops.
                break
            f.add_region(region)

    @staticmethod
    def add_unconformity_above(model: GeologicalModel, feature: GeologicalFeature) -> None:
        """Add a region that excludes values above the unconformity.

        This prevents the value from being
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
    def add_unconformity(
        model: GeologicalModel, feature: GeologicalFeature, value: float, index: int | None = None
    ) -> UnconformityFeature:
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
            if f.type == FeatureType.UNCONFORMITY or f.type == FeatureType.DOMAINFAULT:
                # A domain fault is itself a boundary that already separates
                # everything below it from what's built afterwards (see
                # add_domain_fault_below/above) -- so it's a stopping point
                # here too, the same as an unconformity, rather than
                # something this unrelated unconformity should crop.
                logger.debug(f"Reached unconformity/domain fault {f.name}")
                break
            logger.debug(f"Adding {uc_feature.name} as unconformity to {f.name}")
            if f.type == FeatureType.FAULT or f.type == FeatureType.INACTIVEFAULT:
                continue
            if f == feature:
                continue
            else:
                f.add_region(uc_feature)
        # the faults that are older than the unconformity end on it
        for f in list(model.features):
            if f.type not in (FeatureType.FAULT, FeatureType.INACTIVEFAULT):
                continue
            if FeatureRelationshipManager.fault_is_older_than(model, f, feature):
                FeatureRelationshipManager.truncate_fault(model, f, uc_feature)
        # now add the unconformity to the feature list
        model._add_feature(uc_feature, index=index)
        return uc_feature

    @staticmethod
    def _displaced_by(feature: GeologicalFeature) -> list:
        """The faults that displace a feature (read from its builder)."""
        builder = getattr(feature, "builder", None)
        builders = getattr(builder, "builders", None)
        if builders is not None:
            # a structural frame: the three coordinates have the same faults
            builder = builders[0]
        faults = getattr(builder, "faults", None)
        if faults is None:
            faults = getattr(feature, "faults", [])
        return list(faults)

    @staticmethod
    def fault_is_older_than(
        model: GeologicalModel, fault: GeologicalFeature, feature: GeologicalFeature
    ) -> bool:
        """Whether a fault is older than the unconformity at the base of a feature.

        Call it before the unconformity is added to the model.

        The features of a model are in the order in which they were added,
        and a new feature is displaced by every fault that the model has at
        that moment, unless its ``faults`` argument gives another list. An
        unconformity truncates the features that the model has back to the
        last unconformity or domain fault; the features before that one lie
        above an earlier unconformity and are younger. So the rule is:

        A fault is **younger** than the unconformity if it displaces
        ``feature`` (the series whose isovalue is the unconformity), or a
        feature that lies above an earlier unconformity. Every other fault
        is **older**: a fault added after ``feature``, and a fault added
        before it that the ``faults`` argument of ``feature`` leaves out.

        Parameters
        ----------
        model : GeologicalModel
        fault : FaultSegment
        feature : GeologicalFeature
            the series that defines the unconformity

        Returns
        -------
        bool
        """
        displaced_by = FeatureRelationshipManager._displaced_by
        if any(d is fault for d in displaced_by(feature)):
            return False
        # the features that this unconformity does not truncate: all that
        # are before the last unconformity or domain fault
        boundary = None
        for i in reversed(range(len(model.features))):
            if model.features[i].type in (FeatureType.UNCONFORMITY, FeatureType.DOMAINFAULT):
                boundary = i
                break
        if boundary is None:
            return True
        skip = (
            FeatureType.FAULT,
            FeatureType.INACTIVEFAULT,
            FeatureType.UNCONFORMITY,
            FeatureType.ONLAPUNCONFORMITY,
        )
        for f in model.features[:boundary]:
            if f.type not in skip and any(d is fault for d in displaced_by(f)):
                return False
        return True

    @staticmethod
    def truncate_fault(
        model: GeologicalModel, fault: GeologicalFeature, unconformity: UnconformityFeature
    ) -> bool:
        """Limit a fault to the region below an unconformity.

        The fault frame gets the region of the unconformity: the frame, the
        fault surface and the displacement exist only where the series of
        the unconformity is below its value. A point above the unconformity
        is not moved by the fault, so the fault does not displace the
        features above it. The features that the fault displaces are marked
        to be built again.

        Returns
        -------
        bool
            False if the fault already has this region.
        """
        for r in fault.regions:
            if (
                isinstance(r, UnconformityFeature)
                and r.parent is unconformity.parent
                and r.value == unconformity.value
                and r.sign == unconformity.sign
            ):
                return False
        logger.info(f"Adding {unconformity.name} as unconformity to fault {fault.name}")
        fault.add_region(unconformity)
        for f in model.features:
            if any(d is fault for d in FeatureRelationshipManager._displaced_by(f)):
                builder = getattr(f, "builder", None)
                if builder is not None:
                    builder.set_not_up_to_date(fault)
        return True

    @staticmethod
    def add_unconformities_above_fault(
        model: GeologicalModel, fault: GeologicalFeature, index: int | None = None
    ) -> None:
        """Truncate a new fault by the unconformities that are younger than it.

        A fault that is added after an unconformity is older than it. The
        fault takes the region of every erosional unconformity that is
        before its place in the model (``index``, or the end of the list).
        Onlap unconformities are not used.
        """
        features = model.features if index is None else model.features[:index]
        for f in reversed(features):
            if f.type == FeatureType.UNCONFORMITY:
                FeatureRelationshipManager.truncate_fault(model, fault, f)

    @staticmethod
    def add_onlap_unconformity(
        model: GeologicalModel, feature: GeologicalFeature, value: float, index: int | None = None
    ) -> UnconformityFeature | None:
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
        uc_feature = UnconformityFeature(feature, value, False, onlap=True)
        feature.add_region(uc_feature.inverse())
        for f in reversed(model.features):
            if f.type in (
                FeatureType.UNCONFORMITY,
                FeatureType.ONLAPUNCONFORMITY,
                FeatureType.DOMAINFAULT,
            ):
                logger.debug(f"Reached unconformity {f.name}")
                break
            if f.type == FeatureType.FAULT or f.type == FeatureType.INACTIVEFAULT:
                continue
            if f != feature:
                f.add_region(uc_feature)
        model._add_feature(uc_feature.inverse(), index=index)

        return uc_feature

    @staticmethod
    def relations(model: GeologicalModel) -> dict:
        """Read the relations between the features of a model.

        The relations are read from the builders, which hold them from the
        moment a feature is created, so the result does not depend on whether
        the model has been built. Nothing is changed.

        Parameters
        ----------
        model : GeologicalModel

        Returns
        -------
        dict
            ``events``: the features in the order they were added
            (``{"name", "type"}``, ``type`` a :class:`FeatureType` name). A
            feature added later is older; an unconformity entry follows the
            series that defines it.

            ``faults``: for each fault, ``displaced_by`` (the faults that
            displace its frame, youngest first), ``abuts`` (``{"fault", "side"}``:
            side +1 if the fault exists where the host gx > 0, -1 where it is
            < 0), ``truncated_by`` (see below) and ``displacement``.

            ``features``: for each other feature that is not an unconformity,
            ``type``, ``faults`` (the faults that displace it, in restore
            order: youngest first) and ``truncated_by``.

            ``unconformities``: for each unconformity, ``feature`` (the series
            whose isovalue it is), ``value``, ``onlap`` and ``truncates`` (the
            features and faults whose ``truncated_by`` names it).

            A ``truncated_by`` entry is ``{"unconformity", "feature", "value",
            "keep"}``: the feature exists only where ``feature`` is ``<= value``
            (``keep == "below"``) or ``>= value`` (``"above"``).
        """
        from ...modelling.features.fault import FaultSegment
        from ...utils.regions import NegativeRegion, PositiveRegion

        def faults_of(feature):
            builder = getattr(feature, "builder", None)
            faults = getattr(builder, "faults", None)
            if faults is None:
                faults = feature.faults
            return [f.name for f in faults]

        def truncations(regions):
            out = []
            for r in regions:
                if isinstance(r, UnconformityFeature):
                    out.append(
                        {
                            "unconformity": r.name,
                            "feature": r.parent.name,
                            "value": float(r.value),
                            "keep": "below" if r.sign else "above",
                        }
                    )
            return out

        events, faults, features, unconformities = [], {}, {}, {}
        for feature in model.features:
            events.append({"name": feature.name, "type": FeatureType(feature.type).name})
            if isinstance(feature, FaultSegment):
                abuts = []
                for host, region in feature.abut.items():
                    if isinstance(region, PositiveRegion):
                        side = 1
                    elif isinstance(region, NegativeRegion):
                        side = -1
                    else:
                        side = 0
                    abuts.append({"fault": host, "side": side})
                faults[feature.name] = {
                    "displaced_by": faults_of(feature[0]),
                    "abuts": abuts,
                    "truncated_by": truncations(feature[0].regions),
                    "displacement": float(feature.displacement),
                }
            elif isinstance(feature, UnconformityFeature):
                unconformities[feature.name] = {
                    "feature": feature.parent.name,
                    "value": float(feature.value),
                    "onlap": feature.type == FeatureType.ONLAPUNCONFORMITY,
                    "truncates": [],
                }
            else:
                features[feature.name] = {
                    "type": FeatureType(feature.type).name,
                    "faults": faults_of(feature),
                    "truncated_by": truncations(feature.regions),
                }
        for name, entry in list(faults.items()) + list(features.items()):
            for t in entry["truncated_by"]:
                for u in unconformities.values():
                    if u["feature"] == t["feature"] and u["value"] == t["value"]:
                        u["truncates"].append(name)
        return {
            "events": events,
            "faults": faults,
            "features": features,
            "unconformities": unconformities,
        }
