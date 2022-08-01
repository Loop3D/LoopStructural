from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.modelling.features import FeatureType

import numpy as np


class UnconformityFeature(GeologicalFeature):
    """ """

    def __init__(self, feature: GeologicalFeature, value: float, sign=True):
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
        self.sign = sign

    def inverse(self):
        uc = UnconformityFeature(self, self.value, sign=not self.sign)
        uc.name = self.name + "_inverse"
        return uc

    def evaluate(self, pos: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        pos : numpy array
            locations to evaluate whether below or above unconformity

        Returns
        -------
        np.ndarray.dtype(bool)
            true if above the unconformity, false if below
        """
        if self.sign:
            return self.evaluate_value(pos) < self.value + self.model.epsilon
        if not self.sign:
            return self.evaluate_value(pos) > self.value - self.model.epsilon

    def __call__(self, pos) -> np.ndarray:
        return self.evaluate(pos)
