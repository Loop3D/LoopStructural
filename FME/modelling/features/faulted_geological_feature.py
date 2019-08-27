import numpy as np

from FME.modelling.features.geological_feature import GeologicalFeature
from FME.modelling.scalar_field import ScalarField


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
        # create a geological feature where the property is evaluated on the support nodes
        super().__init__(self.parent_feature.name + "_faulted", ScalarField.from_node_values(
             self.parent_feature.support.support, 
            self.parent_feature.name+'_faulted', 
            self.evaluate_value(self.parent_feature.support.support.nodes)))
    def update(self):
        self.parent_feature.update()
        self.support.update_property(self.evaluate_value(self.parent_feature.support.support.nodes))

    def evaluate_value(self, locations):
        """
        calculate the value of the geological feature at the xyz
        """
        locations = self.fault.apply_to_points(locations)
        return self.parent_feature.evaluate_value(locations)


    def evaluate_gradient(self, locations):
        locations = self.fault.apply_to_points(locations)
        return self.parent_feature.evaluate_gradient(locations)

    def mean(self):
        return np.nanmean(self.support.get_node_values())

    def min(self):
        return np.nanmin(self.support.get_node_values())

    def max(self):
        return np.nanmax(self.support.get_node_values())