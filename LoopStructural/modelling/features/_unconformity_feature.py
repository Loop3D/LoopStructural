from LoopStructural.modelling.features import BaseFeature


class UnconformityFeature(BaseFeature):
    """ """

    def __init__(self, feature, value):
        """

        Parameters
        ----------
        feature
        value
        """
        BaseFeature.__init__(
            self,
            f"{feature.name}_unconformity",
            feature.model,
            feature.faults,
            feature.regions,
            feature.builder,
        )
        self.feature = feature
        self.value = value

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
        return self.feature.evaluate_value(pos) < self.value

    def evaluate_value(self, pos):
        """

        Parameters
        ----------
        pos : numpy array
            locations to evaluate the value of the base geological feature

        Returns
        -------

        """
        return self.feature.evaluate_value(pos)

    def evaluate_gradient(self, pos):
        """

        Parameters
        ----------
        pos : numpy array
            location to evaluate the gradient of the base geological feature

        Returns
        -------

        """
        return self.feature.evaluate_gradient(pos)

    def min(self):
        return self.feature.min()

    def max(self):
        return self.feature.max()
