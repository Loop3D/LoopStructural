class GeologicalFeature:
    def __init__(self):
        self.data = []
        self.interpolator = None
        self.scalar_field = None
        self.associations = {}

    def add_data(self,data):
        """

        :param data:
        :return:
        """
        self.data.append(data)

    def link_feature(self,feature,association):
        self.assocations[association] = feature
        feature.link_feature(self,association)

    def attach_interpolator(self, interpolator):
        self.interpolator = interpolator

    def link_scalar_field(self,scalar_field):
        self.scalar_field = scalar_field

    def evaluate(self,evaluation_points):
        self.scalar_field.evaluate_value(evaluation_points)

    def evaluate_gradient(self,locations):
        self.support.evaluate_gradient(locations)

    def update(self):

        pass


class FaultedGeologicalFeature(GeologicalFeature):
    def __init__(self,feature,fault):
        self.parent_feature = feature
        hw_p, fw_p  = fault.apply_fault_to_support(self.support)
        fault.appy_fault_to_data(data)
        #run interpolator
        self.parent_feature.update()
        hw_v = self.parent_feature.evaluate(hw_p)
        fw_v = self.parent_feature.evaluate(fw_v)

        self.hw_feature = GeologicalFeature()
        self.fw_feature = GeologicalFeature()

        self.hw_feature.link_support()
        #create a new geological feature from



        self.parent_feature.d
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

