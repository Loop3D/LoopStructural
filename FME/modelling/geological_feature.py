from .scalar_field import TetrahedralMeshScalarField


class GeologicalFeature:
    """
    Geological feature is a
    """
    def __init__(self,**kwargs):
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

    def link_feature(self, feature, association):
        self.assocations[association] = feature
        feature.link_feature(self, association)

    def attach_interpolator(self, interpolator):
        self.interpolator = interpolator

    def link_scalar_field(self, scalar_field):
        self.scalar_field = scalar_field

    def evaluate_value(self, evaluation_points):
        self.scalar_field.evaluate_value(evaluation_points)

    def evaluate_gradient(self, locations):
        self.scalar_field.evaluate_gradient(locations)

    def update(self):
        if self.interpolator is None:
            return
        self.interpolator.run()
        self.scalar_field
        pass


class FaultedGeologicalFeature(GeologicalFeature):
    """
    Creates a geological feature that has been faulted using the geological feature representing the fault surface
    and another feature representing the surface pre faulting.
    """
    def __init__(self, feature, fault):
        self.parent_feature = feature
        hw_p, fw_p  = fault.apply_fault_to_support(self.support)
        fault.appy_fault_to_data(data)
        # run interpolator
        self.parent_feature.update()
        hw_v = self.parent_feature.evaluate(hw_p)
        fw_v = self.parent_feature.evaluate(fw_v)

        self.hw_feature = GeologicalFeature()
        self.fw_feature = GeologicalFeature()
        hw_sf = TetrahedralMeshScalarField(self.parent_feature.scalar_field.support, nodevalues=hw_v)
        fw_sf = TetrahedralMeshScalarField(self.parent_feature.scalar_field.support, nodevalues=fw_v)

        self.hw_feature.link_scalar_field(hw_sf)
        self.fw_feature.link_scalar_field(fw_sf)
        # create a new geological feature from

        self.fault = fault

    def evaluate_value(self, locations):
        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros(locations.shape[0])
        evaluated[hangingwall] = self.hangingwall.evaluate(locations[hangingwall])
        evaluated[footwall] = self.footwall.evaluate(locations[footwall])
        return evaluated

    def evaluate_gradient(self, locations):
        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros((locations.shape[0],3))
        evaluated[hangingwall] = self.hangingwall.evaluate(locations[hangingwall])
        evaluated[footwall] = self.footwall.evaluate(locations[footwall])
        return evaluated

