import numpy as np
from LoopStructural.modelling.core.geological_points import IPoint, GPoint, TPoint

import logging
logger = logging.getLogger(__name__)


class GeologicalInterpolator:

    def __init__(self):
        """
        This class is the base class for a geological interpolator and contains all of the
        main interface functions. Any class that is inheriting from this should be callable
        by using any of these functions. This will enable interpolators to be interchanged.
        """
        self.p_i = [] #interface points #   TODO create data container
        self.p_g = [] #gradient points
        self.p_t = [] #tangent points
        self.p_n = [] #norm constraints
        self.n_i = 0
        self.n_g = 0
        self.n_t = 0
        self.n_n = 0
        self.type = 'undefined'
        self.up_to_date = False
        self.constraints = []
        self.propertyname = 'defaultproperty'

    def set_property_name(self,name):
        """
        Set the name of the interpolated property
        Parameters
        ----------
        name string
            name

        Returns
        -------

        """
        self.propertyname = name

    def add_strike_dip_and_value(self, pos, strike, dip, val):
        """
        Add a gradient and value constraint at a location gradient is in the form of strike and dip with the rh thumb

        Parameters
        ----------
        pos
        strike
        dip
        val

        Returns
        -------

        """
        self.n_g +=1
        self.p_g.append(GPoint(pos,strike,dip))
        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))
        self.up_to_date = False

    def add_point(self,pos,val):
        """
        Add a value constraint to the interpolator
        Parameters
        ----------
        pos
        val

        Returns
        -------

        """

        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))

    def add_planar_constraint(self, position, vector, norm=True):
        """
        Add a gradient constraint to the interpolator where the gradient is defined by a normal vector

        Parameters
        ----------
        pos
        val
        norm bool
            whether to constrain the magnitude and directon or only direction

        Returns
        -------

        """
        if norm:
            self.n_n = self.n_n+1
            self.p_n.append(GPoint(position,vector))
        elif norm is False:
            self.n_g = self.n_g+1
            self.p_g.append(GPoint(position,vector))
        self.up_to_date = False

    def add_strike_and_dip(self,pos,s,d, norm=True):
        """
        Add gradient constraint to the interpolator where the gradient is defined by strike and dip

        Parameters
        ----------
        pos
        s
        d
        norm

        Returns
        -------

        """
        self.n_g += 1
        self.p_g.append(GPoint(pos,s,d))
        self.up_to_date = False

    def add_tangent_constraint(self,pos,val):
        """
        Add tangent constraint to the interpolator where the tangent is described by a vector
        :param pos:
        :param val:
        :return:
        """
        self.n_t = self.n_t + 1
        self.p_t.append(TPoint(pos,val))
        self.up_to_date = False

    def add_tangent_constraint_angle(self,pos,s,d):
        """
        Add tangent constraint to the interpolator where the trangent is described by the strike and dip
        :param pos:
        :param s:
        :param d:
        :return:
        """
        self.n_t = self.n_t + 1
        self.p_t.append(TPoint(pos,s,d))
        self.up_to_date = False

    def add_data(self, data):
        """
        Adds a GeologicalData object to the interpolator
        :param data:
        :return:
        """
        if data.type == 'GPoint':
            if not data.norm:
                self.p_g.append(data)
                self.n_g+=1
            if data.norm:
                self.p_n.append(data)
                self.n_n+=1
            self.up_to_date = False
            return
        if data.type == 'IPoint':
            self.p_i.append(data)
            self.n_i+=1
            self.up_to_date = False
            return
        if data.type == 'TPoint':
            self.p_t.append(data)
            self.n_t+=1
            self.up_to_date = False
            return
        else:
            print("Did not add data", data.type)

    def get_value_constraints(self):
        """
        Getter for all active control points
        :return: numpy array Nx4 where 0:3 are the position
        """
        points = np.zeros((self.n_i,4))#array
        for i in range(self.n_i):
            points[i,:3] = self.p_i[i].pos
            points[i,3] = self.p_i[i].val
        return points

    def get_gradient_constraints(self):
        """
        Getter for all gradient control points
        :return: numpy array Nx6 where 0:3 are pos and 3:5 are vector
        """
        points = np.zeros((self.n_g,6))  # array
        for i in range(self.n_g):
            points[i,:3] = self.p_g[i].pos
            points[i,3:] = self.p_g[i].vec
        return points

    def get_tangent_constraints(self):
        points = np.zeros((self.n_t,6))  # array
        for i in range(self.n_t):
            points[i,:3] = self.p_t[i].pos
            points[i,3:] = self.p_t[i].dir
        return points

    def get_norm_constraints(self):
        points = np.zeros((self.n_n,6))
        for i in range(self.n_n):
            points[i,:3] = self.p_n[i].pos
            points[i,3:] = self.p_n[i].vec
        return points

    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        self._setup_interpolator(**kwargs)

    def solve_system(self,**kwargs):
        """
        Solves the interpolation equations
        """
        self._solve(**kwargs)
        self.up_to_date = True

