import numpy as np

from FME.modelling.features.geological_feature import GeologicalFeature
from FME.modelling.scalar_field import TetrahedralMeshScalarField


class FaultedGeologicalFeature(GeologicalFeature):
    """
    Creates a geological feature that has been faulted using the geological feature representing the fault
    surface and another feature representing the surface pre faulting.
    """
    def __init__(self, feature, fault):
        super().__init__(feature.name+"_faulted", feature.support)
        self.parent_feature = feature
        hw_p, fw_p, hw_m, fw_m = fault.apply_fault_to_support(self.support)
        # fault.appy_fault_to_data(data)

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
        self.hw_feature = GeologicalFeature(self.name+'_hw', hw_sf)
        self.fw_feature = GeologicalFeature(self.name+'_fw', fw_sf)
        self.fault = fault
        #create a continuous geological feature
        faulted_v = np.zeros(self.support.number_of_nodes())
        fault_mask = fault.faultframe.features[0].support.get_node_values()>0
        faulted_v[fault_mask] = hw_v[fault_mask]
        faulted_v[~fault_mask] = fw_v[~fault_mask]
        self.feature = GeologicalFeature(self.name+"_faulted",
                                         TetrahedralMeshScalarField.
                                         from_node_values(self.support.mesh,
                                                          self.name + '_faulted',
                                                          faulted_v))
    def evaluate_value(self, locations):
        """
        calculate the value of the geological feature at the xyz
        """
        hangingwall = self.fault.evaluate(locations)
        footwall = ~self.fault.evaluate(locations)
        evaluated = np.zeros(locations.shape[0])
        evaluated[hangingwall] = self.hw_feature.evaluate_value(locations[hangingwall])
        evaluated[footwall] = self.fw_feature.evaluate_value(locations[footwall])
        return evaluated

    def evaluate_gradient(self, locations):

        hangingwall = self.fault.evaluate(locations) > 0
        footwall = self.fault.evaluate(locations) < 0
        evaluated = np.zeros((locations.shape[0], 3))
        evaluated[hangingwall] = self.hw_feature.evaluate_value(locations[hangingwall])
        evaluated[footwall] = self.fw_feature.evaluate_value(locations[footwall])
        return evaluated

    def mean(self):
        return np.nanmean(self.support.get_node_values())

    def min(self):
        return np.nanmin(self.support.get_node_values())

    def max(self):
        return np.nanmax(self.support.get_node_values())