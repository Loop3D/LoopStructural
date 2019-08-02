from .scalar_field import TetrahedralMeshScalarField
from ..modelling.geological_points import GPoint, IPoint, TPoint
import numpy as np


class GeologicalFeatureBuilder:
    def __init__(self, interpolator, **kwargs):
        self.interpolator = interpolator
        if 'name' in kwargs:
            self.name = kwargs['name']
        if 'region' in kwargs:
            self.region = kwargs['region']
        self.data = []

    def add_strike_dip_and_value(self, pos, strike, dip, val):
        self.data.append(GPoint(pos, strike, dip))
        self.data.append(IPoint(pos, val))

    def add_point(self, pos, val):
        self.data.append(IPoint(pos, val))

    def add_planar_constraint(self, pos, val):
        self.data.append(GPoint(pos, val))

    def add_strike_and_dip(self, pos, s, d):
        self.data.append(GPoint(pos, s, d))

    def add_tangent_constraint(self, pos, val):
        self.data.append(TPoint(pos, val))

    def add_tangent_constraint_angle(self, pos, s, d):
        self.data.append(TPoint(pos, s, d))

    def build(self, solver='cg', **kwargs):
        for d in self.data:
            self.interpolator.add_data(d)
        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(solver=solver)
        return GeologicalFeature(self.name,
                                 TetrahedralMeshScalarField.from_interpolator(self.interpolator))


class GeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    modle. For example foliations, fault planes, fold rotation angles etc. The feature has a support
    which 
    """
    def __init__(self, name, support):
        """

        :param age:
        :param name:
        :param support:
        """
        self.name = name
        self.support = support
        self.ndim = 1
        self.data = []

    def evaluate_value(self, evaluation_points):
        return self.support.evaluate_value(evaluation_points)

    def evaluate_gradient(self, locations):
        return self.support.evaluate_gradient(locations)

    def mean_property_value(self):
        return np.nanmean(self.support.get_node_values())

    def min_property_value(self):
        return np.nanmin(self.support.get_node_values())

    def max_property_value(self):
        return np.nanmax(self.support.get_node_values())


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

    def mean_property_value(self):
        return np.nanmean(self.support.get_node_values())

    def min_property_value(self):
        return np.nanmin(self.support.get_node_values())

    def max_property_value(self):
        return np.nanmax(self.support.get_node_values())


class CompositeGeologicalFeature(GeologicalFeature):
    def __init__(self, geological_feature_a, region_a, geological_feature_b, region_b):
        self.geological_feature_a = geological_feature_a
        self.region_a = region_a

        self.geological_feature_b = geological_feature_b
        self.region_b = region_b

    def evaluate_gradient(self, locations):
        
        pass

    def evaluate_value(self, evaluation_points):

        pass

    def mean_property_value(self):
        return np.nanmean(self.support.get_node_values())

    def min_property_value(self):
        return np.nanmin(self.support.get_node_values())

    def max_property_value(self):
        return np.nanmax(self.support.get_node_values())


class CrossProductGeologicalFeature(GeologicalFeature):
    def __init__(self, name, geological_feature_a, geological_feature_b):
        super().__init__(name+"_faulted", geological_feature_a.support)
        self.geological_feature_a = geological_feature_a
        self.geological_feature_b = geological_feature_b
    def evaluate_gradient(self, locations):
        return np.cross(self.geological_feature_a.evaluate_gradient(locations),
                        self.geological_feature_b.evaluate_gradient(locations),
                        axisa=1,
                        axisb=1)
    def evaluate_value(self, evaluation_points):
        return 0.

    def mean_property_value(self):
        return 0.

    def min_property_value(self):
        return 0.

    def max_property_value(self):
        return 0.
