class RegionFeature:
    def __init__(self, function):
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
        return 'region'