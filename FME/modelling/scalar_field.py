class ScalarField:
    def __init__(self,support,property_name):
        self.support = support
        self.property_name = property_name

    def evaluate_value(self,evaluation_points):
        return self.support.evaluate_property(self.property_name,evaluation_points)

    def evaluate_gradient(self,evaluation_points):
        return self.support.evaluate_gradient(self.property_name,evaluation_points)