from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.modelling.features import FeatureType


class UnconformityFeature(GeologicalFeature):
    """ """

    def __init__(self, feature: GeologicalFeature, value: float):
        """

        Parameters
        ----------
        feature
        value
        """
        # create a shallow(ish) copy of the geological feature
        # just don't link the regions
        GeologicalFeature.__init__(
            self,
            name=f"{feature.name}_unconformity",
            faults=feature.faults,
            regions=[],  # feature.regions.copy(),  # don't want to share regionsbetween unconformity and # feature.regions,
            builder=feature.builder,
            model=feature.model,
            interpolator=feature.interpolator,
        )
        self.value = value
        self.type = FeatureType.UNCONFORMITY

    def evaluate(self, pos):
        """

        Parameters
        ----------
        pos : numpy array
            locations to evaluate whether below or above unconformity

        Returns
        -------
        boolean
            true if above the unconformity, false if below
        """
        return self.evaluate_value(pos) < self.value
