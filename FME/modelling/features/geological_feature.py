from FME.modelling.scalar_field import ScalarField
from FME.modelling.geological_points import GPoint, IPoint, TPoint
import numpy as np


class GeologicalFeatureInterpolator:
    """
    A builder for a GeologicalFeature
    """
    def __init__(self, interpolator, **kwargs):
        self.interpolator = interpolator
        self.name = "UnnamedFeature"
        if 'name' in kwargs:
            self.name = kwargs['name']
            self.interpolator.set_property_name(self.name)
        if 'region' in kwargs:
            self.region = kwargs['region']
        self.data = []

    def update(self):
        pass

    def add_strike_dip_and_value(self, pos, strike, dip, val):
        """

        Parameters
        ----------
        pos
        strike
        dip
        val

        Returns
        -------

        """
        self.data.append(GPoint(pos, strike, dip))
        self.data.append(IPoint(pos, val))

    def add_point(self, pos, val):
        """

        Parameters
        ----------
        pos
        val

        Returns
        -------

        """
        self.data.append(IPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_planar_constraint(self, pos, val):
        """

        Parameters
        ----------
        pos
        val

        Returns
        -------

        """
        self.data.append(GPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_strike_and_dip(self, pos, s, d):
        """

        Parameters
        ----------
        pos
        s
        d

        Returns
        -------

        """
        self.data.append(GPoint.from_strike_and_dip(pos, s, d))
        self.interpolator.add_data(self.data[-1])

    def add_plunge_and_plunge_dir(self,pos,plunge,plunge_dir):
        """

        Parameters
        ----------
        pos
        plunge
        plunge_dir

        Returns
        -------

        """
        self.data.append(GPoint.from_plunge_plunge_dir(pos,plunge,plunge_dir))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint(self, pos, val):
        """

        Parameters
        ----------
        pos
        val

        Returns
        -------

        """
        self.data.append(TPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint_angle(self, pos, s, d):
        """

        Parameters
        ----------
        pos
        s
        d

        Returns
        -------

        """
        self.data.append(TPoint(pos, s, d))
        self.interpolator.add_data(self.data[-1])

    def build(self, solver='cg', **kwargs):
        """

        Parameters
        ----------
        solver
        kwargs

        Returns
        -------

        """
        # for d in self.data:
        #     self.interpolator.add_data(d)
        # we can add a fold to the interpolator if the interpolator is a fold interpolator
        # pass the dict with weights as kwargs to the fold interpolator
        if "fold" in kwargs and "fold_weights" in kwargs:
            self.interpolator.update_fold(kwargs['fold'])
            self.interpolator.add_fold_constraints(**kwargs['fold_weights'])
            kwargs['cg'] = False
        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(solver=solver)
        return GeologicalFeature(self.name,
                                 ScalarField.from_interpolator(self.interpolator),
                                 builder=self, data=self.data)


class GeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    modle. For example foliations, fault planes, fold rotation angles etc. The feature has a support
    which 
    """
    def __init__(self, name, support, builder = None, data = None):
        """

        :param age:
        :param name:
        :param support:
        """
        self.name = name
        self.support = support
        self.ndim = 1
        self.data = data
        self.builder = builder

    def set_builder(self, builder):
        """

        Parameters
        ----------
        builder

        Returns
        -------

        """
        self.builder = builder

    def evaluate_value(self, evaluation_points):
        """

        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        return self.support.evaluate_value(evaluation_points)

    def evaluate_gradient(self, locations):
        """

        Parameters
        ----------
        locations

        Returns
        -------

        """
        return self.support.evaluate_gradient(locations)

    def mean(self):
        """

        Returns
        -------

        """
        return np.nanmean(self.support.get_node_values())

    def min(self):
        """

        Returns
        -------

        """
        return np.nanmin(self.support.get_node_values())

    def max(self):
        """

        Returns
        -------

        """
        return np.nanmax(self.support.get_node_values())

    def update(self):
        """

        Returns
        -------

        """
        self.support.interpolator.up_to_date = False
        self.support.interpolator.update()

