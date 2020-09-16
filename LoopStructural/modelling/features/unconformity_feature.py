"""

"""
class UnconformityFeature:
    """
    
    """
    def __init__(self, feature, value):
        """

        Parameters
        ----------
        feature
        value
        """
        self.feature = feature
        self.value = value
        self.type = 'unconformity'
        self.name = '{}_unconformity'.format(feature.name)
        self.builder = self.feature.builder
    def add_region(self, region):
        # self.feature.add_region(region)
        pass
    def set_model(self, model):
        self.model = model

    def evaluate(self,pos):
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

    def evaluate_value(self,pos):
        """

        Parameters
        ----------
        pos : numpy array
            locations to evaluate the value of the base geological feature

        Returns
        -------

        """
        return self.feature.evaluate_value(pos)

    def evaluate_gradient(self,pos):
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