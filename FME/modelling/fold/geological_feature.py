class GeologicalFeature:
    def __init__(self):
        self.data = []
        self.interpolator = None
        self.support = None
        self.geological_structures = {}
    def add_data(self,data):
        self.data.append(data)
    def link_interpolator(self,interpolator):
        self.interpolator = interpolator
    def link_support(self,support):
        self.support = support
    def evaluate(self,locations):
        self.support.evaluate(locations)
    def evaluate_gradient(self,locations):
        self.support.evaluate_gradient(locations)
class FaultedGeologicalFeature(GeologicalFeature):
    def __init__(self,feature,fault):
        self.parent_feature = feature
        self.hangingwall = hangingwall
        self.footwall = footwall
        self.fault = fault
    def evaluate(self,locations):
        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros(locations.shape[0])
        evaluated[hangingwall] = self.hangingwall.evaluate(locations[hangingwall])
        evaluated[footwall] = self.footwall.evaluate(locations[footwall])

        return evaluated

