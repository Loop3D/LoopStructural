from FME.modelling.scalar_field import TetrahedralMeshScalarField
from FME.modelling.geological_points import GPoint, IPoint, TPoint
import numpy as np


class GeologicalFeatureInterpolator:
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
        self.interpolator.add_data(self.data[-1])

    def add_planar_constraint(self, pos, val):
        self.data.append(GPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_strike_and_dip(self, pos, s, d):
        self.data.append(GPoint(pos, s, d))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint(self, pos, val):
        self.data.append(TPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint_angle(self, pos, s, d):
        self.data.append(TPoint(pos, s, d))
        self.interpolator.add_data(self.data[-1])

    def build(self, solver='cg', **kwargs):
        # for d in self.data:
        #     self.interpolator.add_data(d)
        # we can add a fold to the interpolator if the interpolator is a fold interpolator
        # pass the dict with weights as kwargs to the fold interpolator
        if "fold" in kwargs and "fold_weights" in kwargs:
            self.interpolator.update_fold(kwargs['fold'])
            self.interpolator.add_fold_constraints(**kwargs['fold_weights'])

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

    def mean(self):
        return np.nanmean(self.support.get_node_values())

    def min(self):
        return np.nanmin(self.support.get_node_values())

    def max(self):
        return np.nanmax(self.support.get_node_values())


