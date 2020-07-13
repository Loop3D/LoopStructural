class RegionFeature:
    """
    """
    def __init__(self, function):
        """
        Create a GeologicalFeature to represent a region in a model
        The region is defined by a boolean function on position.

        Parameters
        ----------
        function lambda function
            lambda function true inside region, false outside region
        
        """
        self.function = function
        self.name = 'region'

    def evaluate_value(self, pos):
        return self.function(pos).astype(float)

    def mean(self):
        return 0

    def max(self):
        return 1

    def min(self):
        return -1

    def name(self):
        return self.name

# class VectorFeature:
#     def __init__(self, function)
