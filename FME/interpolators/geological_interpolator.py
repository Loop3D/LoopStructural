import numpy as np
from FME.modelling.geological_points import IPoint, GPoint, TPoint


class GeologicalInterpolator:
    """
    This class is the base class for a geological interpolator and contains all of the
    main interface functions. Any class that is inheriting from this should be callable
    by using any of these functions. This will enable interpolators to be interchanged.
    """
    def __init__(self,**kwargs):
        """
        Default constructor requires no arguments
        :param kwargs: 'itype' what type of geological feature is being interpolated
        """
        self.p_i = [] #interface points #TODO create data container
        self.p_g = [] #gradeint points
        self.p_t = [] #tangent points
        self.n_i = 0
        self.n_g = 0
        self.n_t = 0
        self.type = 'undefined'
        if 'itype' in kwargs:
            self.type = kwargs['itype']
        self.up_to_date = False
        self.constraints = []
        self.headings = ["Constraint Type","Number of constraints", "Per Constraint Weighting"]
    def add_strike_dip_and_value(self,pos,strike,dip,val):
        """
        Add a gradient and value constraint at a location gradient is in the form of strike and dip with the rh thumb
        rule
        :param pos:
        :param strike:
        :param dip:
        :param val:
        :return:
        """
        self.n_g +=1
        self.p_g.append(GPoint(pos,strike,dip))
        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))
        self.up_to_date = False

    def add_point(self,pos,val):
        """
        Add interface point to the interpolator
        :param pos:
        :param val:
        :return:
        """
        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))

    def add_planar_constraint(self,pos,val):
        """
        Add a gradient constraint to the interpolator where the gradient is defined by a normal vector
        """
        self.n_g = self.n_g+1
        self.p_g.append(GPoint(pos,val))
        self.up_to_date = False


    def add_strike_and_dip(self,pos,s,d):
        """
        Add gradient constraint to the interpolator where the gradient is defined by strike and dip
        :param pos:
        :param s:
        :param d:
        :return:
        """
        self.n_g +=1
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

    def add_data(self,data):
        """
        Adds a GeologicalData object to the interpolator
        :param data:
        :return:
        """
        if type(data) == GPoint:
            self.p_g.append(data)
            self.n_g+=1
        if type(data) == IPoint:
            self.p_i.append(data)
            self.n_i+=1
        if type(data) == TPoint:
            self.p_t.append(data)
            self.n_t+=1
        self.up_to_date = False

    def get_control_points(self):
        """
        Getter for all active control points
        :return: numpy array Nx4 where 0:3 are the position
        """
        points = np.zeros((self.n_i,4))#array
        for i in range(self.n_i):
            points[i,:3] = self.p_i[i].pos
            points[i,3] = self.p_i[i].val
        return points

    def get_gradient_control(self):
        """
        Getter for all gradient control points
        :return: numpy array Nx6 where 0:3 are pos and 3:5 are vector
        """
        points = np.zeros((self.n_g,6))  # array
        for i in range(self.n_g):
            points[i,:3] = self.p_g[i].pos
            points[i,3:] = self.p_g[i].dir
        return points

    def get_tangent_control(self):
        points = np.zeros((self.n_t,6))  # array
        for i in range(self.n_t):
            points[i,:3] = self.p_t[i].pos
            points[i,3:] = self.p_t[i].dir
        return points

    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        #print(columnar.columnar(self.constraints,self.headings))
        self._setup_interpolator(**kwargs)

    def interpolate_value(self,points):
        """
        Evaluate the interpolator at the points
        :param points:
        :return:
        """
        return self._interpolate_value(points)

    def interpolate_gradient(self,points):
        """
        Evaluate the inteprolator gradient at the points
        :param points:
        :return:
        """
        return self._interpolate_gradient(points)

    def solve_system(self,**kwargs):
        """
        Solves the interpolation equations
        """
        self._solve(**kwargs)
        self.up_to_date = True

