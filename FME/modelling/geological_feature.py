from .scalar_field import TetrahedralMeshScalarField
import numpy as np


class GeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element of a geological
    feature. For example foliations, fault planes, fold rotation angles etc.
    """
    def __init__(self, age, name, support):
        """

        :param age:
        :param name:
        :param support:
        """
        self.age = age
        self.name = name
        self.support = support
        self.ndim = 1
    # @classmethod
    # def from_interpolation_parameters(cls, age, name, support, params):
    #
    #     pass

    def evaluate_value(self, evaluation_points):
        return self.support.evaluate_value(evaluation_points)

    def evaluate_gradient(self, locations):
        return self.support.evaluate_gradient(locations)


class FaultedGeologicalFeature(GeologicalFeature):
    """
    Creates a geological feature that has been faulted using the geological feature representing the fault surface
    and another feature representing the surface pre faulting.
    """
    def __init__(self, feature, fault):
        super().__init__(feature.age, feature.name+"_faulted", feature.support)
        self.parent_feature = feature
        hw_p, fw_p, hw_m, fw_m  = fault.apply_fault_to_support(self.support)
        # fault.appy_fault_to_data(data)
        # run interpolator
        # self.parent_feature.update()
        hw_v = np.zeros(self.support.number_of_nodes())
        fw_v = np.zeros(self.support.number_of_nodes())
        hw_v[:] = np.nan
        fw_v[:] = np.nan
        #
        hw_v[hw_m] = self.parent_feature.evaluate_value(hw_p)
        fw_v[fw_m] = self.parent_feature.evaluate_value(fw_p)
        # print(hw_v)
        # hw_v = self.parent_feature.support.get_node_values()
        hw_sf = TetrahedralMeshScalarField.from_node_values(self.support.mesh, self.name+'_hw', hw_v)
        fw_sf = TetrahedralMeshScalarField.from_node_values(self.support.mesh, self.name+'_fw', fw_v)
        self.hw_feature = GeologicalFeature(self.age, self.name+'_hw', hw_sf)
        self.fw_feature = GeologicalFeature(self.age, self.name+'_fw', fw_sf)
        # create a new geological feature from

        self.fault = fault

    def evaluate_value(self, locations):

        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros(locations.shape[0])
        evaluated[hangingwall] = self.hw_feature.evaluate(locations[hangingwall])
        evaluated[footwall] = self.fw_feature.evaluate(locations[footwall])
        return evaluated

    def evaluate_gradient(self, locations):

        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros((locations.shape[0],3))
        evaluated[hangingwall] = self.hw_feature.evaluate(locations[hangingwall])
        evaluated[footwall] = self.fw_feature.evaluate(locations[footwall])
        return evaluated

