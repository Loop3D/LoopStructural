"""
Base geological interpolator
"""
import logging

import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class GeologicalInterpolator:
    """
    Attributes
    ----------
    data : dict
        a dictionary with np.arrays for gradient, value, normal, tangent data
    """
    def __init__(self):
        """
        This class is the base class for a geological interpolator and contains all of the
        main interface functions. Any class that is inheriting from this should be callable
        by using any of these functions. This will enable interpolators to be interchanged.
        """

        self.data = {'gradient': np.zeros((0,7)),
                     'value' : np.zeros((0,5)),
                     'normal': np.zeros((0,7)),
                     'tangent': np.zeros((0,7)),
                     'interface' : np.zeros((0,5))
                     }
        self.n_g = 0
        self.n_i = 0
        self.n_n = 0
        self.n_t = 0

        self.type = 'undefined'
        self.up_to_date = False
        self.constraints = []
        self.propertyname = 'defaultproperty'
        self.__str = 'Base Geological Interpolator'

    def __str__(self):
        
        return self.__str

    def set_region(self,**kwargs):
        pass

    def set_property_name(self, name):
        """
        Set the name of the interpolated property
        Parameters
        ----------
        name : string
            name of the property to be saved on a mesh

        Returns
        -------

        """
        self.propertyname = name

    def set_value_constraints(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """

        self.data['value'] = points
        self.n_i = points.shape[0]

    def set_gradient_constraints(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        self.n_g = points.shape[0]
        self.data['gradient'] = points

    def set_normal_constraints(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        self.n_n = points.shape[0]
        self.data['normal'] = points

    def set_tangent_constraints(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        self.data['tangent'] = points

    def set_interface_constraints(self, points):
        self.data['interface'] = points

    def get_value_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data['value']

    def get_gradient_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data['gradient']


    def get_tangent_constraints(self):
        """

        Returns
        -------
        numpy array
        """

        return self.data['tangent']

    def get_norm_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data['normal']

    def get_data_locations(self):
        norm = self.get_norm_constraints()
        grad = self.get_gradient_constraints()
        val = self.get_value_constraints()
        return np.vstack([norm[:,:3],grad[:,:3],val[:,:3]])
        
    def get_interface_constraints(self):
        return self.data['interface']
    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        self._setup_interpolator(**kwargs)

    def solve_system(self, **kwargs):
        """
        Solves the interpolation equations
        """
        self._solve(**kwargs)
        self.up_to_date = True

    def update(self):
        return False

    def reset(self):
        """
        Removes all of the data from an interpolator

        Returns
        -------

        """
        self.n_g = 0
        self.n_i = 0
        self.n_n = 0
        self.n_t = 0

