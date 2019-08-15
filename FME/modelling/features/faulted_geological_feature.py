import numpy as np

from FME.modelling.features.geological_feature import GeologicalFeature
from FME.modelling.scalar_field import TetrahedralMeshScalarField


class FaultedGeologicalFeature(GeologicalFeature):
    """
    Creates a geological feature that has been faulted using the geological feature representing the fault
    surface and another feature representing the surface pre faulting.
    """
    def __init__(self, feature, fault):
        # start with the feature to be faulted this is the parent feature
        self.parent_feature = feature
        self.fault = fault
        self.hw_feature = None
        self.fw_feature = None
        self.fault.apply_to_data(self.parent_feature.data)
        self.update_feature()

    def update(self):
        self.parent_feature.update()

    def update_feature(self):
        # determine the hw and fw movements
        hw_p, fw_p, hw_m, fw_m = self.fault.apply_to_support(self.parent_feature.support)
        # TODO this should all be managed by an observer class which links the data
        # to the feature/interpolator and tell the interpolator that it needs to rerun
        # if type(self.parent_feature) == GeologicalFeature:
        self.update()
        # else:
        #     self.parent_feature.update_feature()
        # evaluate the values of the faulted points
        hw_v = np.zeros(self.parent_feature.support.number_of_nodes())
        fw_v = np.zeros(self.parent_feature.support.number_of_nodes())
        hw_v[:] = np.nan
        fw_v[:] = np.nan
        #
        hw_v[hw_m] = self.parent_feature.evaluate_value(hw_p)
        fw_v[fw_m] = self.parent_feature.evaluate_value(fw_p)
        # print(hw_v)
        # hw_v = self.parent_feature.support.get_node_values()
        # create new supports for the hw and fw features
        hw_sf = TetrahedralMeshScalarField.from_node_values(
            self.parent_feature.support.mesh, self.parent_feature.name+'_hw', hw_v)
        fw_sf = TetrahedralMeshScalarField.from_node_values(
            self.parent_feature.support.mesh, self.parent_feature.name+'_fw', fw_v)
        self.hw_feature = GeologicalFeature(self.parent_feature.name+'_hw', hw_sf)
        self.fw_feature = GeologicalFeature(self.parent_feature.name+'_fw', fw_sf)
        # self.fault = fault
        # now initialise the actual feature so it can be called
        #create a continuous geological feature
        faulted_v = np.zeros(self.parent_feature.support.number_of_nodes())
        fault_mask = self.fault.faultframe.features[0].support.get_node_values()>0
        faulted_v[fault_mask] = hw_v[fault_mask]
        faulted_v[~fault_mask] = fw_v[~fault_mask]
        super().__init__(self.parent_feature.name + "_faulted", TetrahedralMeshScalarField.
                          from_node_values(self.parent_feature.support.mesh,
                                           self.parent_feature.name + '_faulted',
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
        hangingwall = self.fault.evaluate(locations)
        footwall = ~self.fault.evaluate(locations)
        evaluated = np.zeros((locations.shape[0], 3))
        evaluated[hangingwall] = self.hw_feature.evaluate_gradient(locations[hangingwall])
        evaluated[footwall] = self.fw_feature.evaluate_gradient(locations[footwall])
        return evaluated

    def mean(self):
        return np.nanmean(self.support.get_node_values())

    def min(self):
        return np.nanmin(self.support.get_node_values())

    def max(self):
        return np.nanmax(self.support.get_node_values())