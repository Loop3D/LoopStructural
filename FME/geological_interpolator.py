import sys
import numpy as np
from .geological_points import IPoint,GPoint,TPoint,IePoint
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
            points[i,3] = self.p_i[i].val
        return points
    def get_gradient_control(self):
        points = np.zeros((self.n_g,6))#array
        for i in range(self.n_g):
            points[i,:3] = self.p_g[i].pos
            points[i,3:] = self.p_g[i].dir
        return points
    def get_tangent_control(self):
        points = np.zeros((self.n_t,6))#array
        for i in range(self.n_t):
            points[i,:3] = self.p_t[i].pos
            points[i,3:] = self.p_t[i].dir
        return points
    def setup_interpolator(self,**kwargs):
        """
        Runs all of the required setting up stuff
        """
        self._setup_interpolator(**kwargs)
    def interpolate_value(self,points):
        return self._interpolate_value(points)
    def interpolate_gradient(self,points):
        return self._interpolate_gradient(points)
    def solve_system(self,**kwargs):
        self._solve(**kwargs)
    def export_data_to_vtk(self,filename):
        from pyevtk.hl import pointsToVTK

        points = self.get_control_points()
        if points.shape[0] > 0:
            x = np.array(points[:,0], copy=True, order='C')
            y = np.array(points[:,1], copy=True, order='C')
            z = np.array(points[:,2], copy=True, order='C')
            v = np.array(points[:,3], copy=True, order='C')
            pointsToVTK(filename+"_value_points",x,y,z,data={"val":v})
        points = self.get_gradient_control()
        if points.shape[0] > 0:
            x = np.array(points[:,0], copy=True, order='C')
            y = np.array(points[:,1], copy=True, order='C')
            z = np.array(points[:,2], copy=True, order='C')
            v = (np.array(points[:,3], copy=True, order='C'),np.array(points[:,4], copy=True, order='C'),np.array(points[:,5], copy=True, order='C'))
            pointsToVTK(filename+"_gradient_points",x,y,z,data={"v":v})
        points = self.get_tangent_control()
        if points.shape[0] > 0:
            x = np.array(points[:,0], copy=True, order='C')
            y = np.array(points[:,1], copy=True, order='C')
            z = np.array(points[:,2], copy=True, order='C')
            v = (np.array(points[:,3], copy=True, order='C'),np.array(points[:,4], copy=True, order='C'),np.array(points[:,5], copy=True, order='C'))#v = np.array(points[:,3:], copy=True, order='C')
           
            pointsToVTK("./"+filename+"_tangent_points",x,y,z,data={"v":v})
