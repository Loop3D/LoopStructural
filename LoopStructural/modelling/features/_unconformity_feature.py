import numpy as np

from ...modelling.features import FeatureType, GeologicalFeature


class UnconformityFeature(GeologicalFeature):
    """ """

    def __init__(self, feature: GeologicalFeature, value: float, sign=True, onlap=False):
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
            name=f"__{feature.name}_unconformity",
            faults=feature.faults,
            regions=[],  # feature.regions.copy(),  # don't want to share regionsbetween unconformity and # feature.regions,
            builder=feature.builder,
            model=feature.model,
        )
        self.value = value
        self.type = FeatureType.UNCONFORMITY if onlap is False else FeatureType.ONLAPUNCONFORMITY
        self.sign = sign
        self.parent = feature

    @property
    def faults(self):
        return self.parent.faults

    def to_json(self):
        json = super().to_json()
        json["value"] = self.value
        json["sign"] = self.sign
        json["parent"] = self.parent.name
        return json

    def inverse(self):
        """Returns an unconformity feature with the sign flipped
        The feature is a shallow copy with the parent being set to
        the parent of this feature

        Returns
        -------
        UnconformityFeature
            a new unconformity feature with the sign of the unconformity reversed
        """
        uc = UnconformityFeature(
            self.parent,
            self.value,
            sign=not self.sign,
            onlap=self.type == FeatureType.ONLAPUNCONFORMITY,
        )
        uc.name = self.name + "_inverse"
        return uc

    # _calculate_mask passes restored_through to a region with this flag
    takes_restored_through = True

    def evaluate(self, pos: np.ndarray, *, block=None, restored_through=()) -> np.ndarray:
        """

        Parameters
        ----------
        pos : numpy array
            locations to evaluate whether below or above unconformity
        block : dict, optional
            fault block in which the unconformity is evaluated
        restored_through : sequence of FaultSegment, optional
            faults whose restoration the points already have: the series of
            the unconformity does not apply them again

        Returns
        -------
        np.ndarray.dtype(bool)
            true if above the unconformity, false if below
        """
        kwargs = {"restored_through": tuple(restored_through)} if restored_through else {}
        if self.sign:
            return self.evaluate_value(pos, block=block, **kwargs) <= self.value
        if not self.sign:
            return self.evaluate_value(pos, block=block, **kwargs) >= self.value

    def __call__(self, pos, *, block=None, restored_through=()) -> np.ndarray:
        return self.evaluate(pos, block=block, restored_through=restored_through)
