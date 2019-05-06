import sys
import numpy as np

from geological_points import IPoint,GPoint,TPoint,IePoint
class GeologicalInterpolator():
    """
    This class is the base class for a geological interpolator and contains all of the
    main interface functions. Any class that is inheriting from this should be callable
    by using any of these functions. This will enable interpolators to be interchanged.
    """
    def __init__(self,**kwargs):
        self.p_i = [] #interface points
        self.p_g = [] #gradeint points
        self.p_t = [] #tangent points
        self.n_i = 0
        self.n_g = 0
        self.n_t = 0
        self.type = 'undefined'
        if 'itype' in kwargs:
            self.type = kwargs['itype']
    def add_strike_dip_and_value(self,pos,strike,dip,val):
        self.n_g +=1
        self.p_g.append(GPoint(pos,strike,dip))
        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))
    def add_point(self,pos,val):
        self.n_i = self.n_i + 1
        self.p_i.append(IPoint(pos,val))
    def add_planar_constraint(self,pos,val):
        self.n_g = self.n_g+1
        self.p_g.append(GPoint(pos,val))
    def add_strike_and_dip(self,pos,s,d):
        self.n_g +=1
        self.p_g.append(GPoint(pos,s,d))
    def add_tangent_constraint(self,pos,val):
        self.n_t = self.n_t + 1
        self.p_t.append(TPoint(pos,val))
    def add_tangent_constraint_angle(self,pos,s,d):
        self.n_t = self.n_t + 1
        self.p_t.append(TPoint(pos,s,d))
    def add_data(self,data):
        if type(data) == GPoint:
            self.p_g.append(data)
            self.n_g+=1
        if type(data) == IPoint:
            self.p_i.append(data)
            self.n_i+=1
        if type(data) == TPoint:
            self.p_t.append(data)
            self.n_t+=1
                
    def fold_event(self,fold):
        self.fold = fold
    def get_control_points(self):
        points = np.zeros((self.n_i,4))#array
        for i in range(self.n_i):
            points[i,:3] = self.p_i[i].pos

