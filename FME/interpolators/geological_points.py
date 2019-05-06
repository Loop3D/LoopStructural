import numpy as np
import math as m
from numpy import linalg as nla
import dist


class Point():
    def __init__(self,pos):
        self.pos = np.array(pos)
        self.orig = np.array(pos)
    def transform(self,t):
        #self.orig = self.pos
        t /= nla.norm(t)
        t[0] = - t[1]
        t[1] = t[0]
        t[2] = t[2]
        #print('t',t)
        self.pos = self.pos+t
    def restore(self):
        self.pos = self.orig    
    def dist(self,point):
        eps = np.finfo('float').eps
        #check points are the same dimensions
        if len(self.pos) != len(point.pos):
            return False
        if self.pos.ndim == 2:# > point.pos.shape[1]:
            pos = np.tile(point.pos,(self.pos.shape[1],1)).T
            return nla.norm(pos - self.pos,axis=1)+eps, pos - self.pos        
        if point.pos.ndim == 2:# > point.pos.shape[1]:
            d = np.zeros((point.pos.shape[1]))
            delta = np.zeros((point.pos.shape[1],3))
            if dist.dist_2d(point.pos.T,self.pos,d,delta):
                #add epsilon to radius to avoid sqrt(0) and /0
                return d+eps, delta.T
            else:
                raise BaseException
        d = 0.0
        d = dist.dist_1d(point.pos,self.pos,d)+eps
        return d,point.pos-self.pos#nla.norm(point.pos - self.pos,axis=0)+eps, point.pos-self.pos        
    def dim(self):
        return len(self.pos)
class GridPoint(Point):
    """
    Point in a grid
    """
    def __init__(self,pos):
        Point.__init__(self,pos)
        self.properties = {}
    def add_property(self,name,value):
        self.properties[name] = value
    def element(self,i):
        self.element=i
class IPoint(Point):
    """
    Interface point
    """
    def __init__(self,pos,val):
        Point.__init__(self,pos)
        self.val = val
    def val(self):
        return self.val
class IePoint(Point):
    """
    Inequality point
    """
    def __init__(self,pos,val,greater=True):
        Point.__init__(self,pos)
        self.val = val
        self.greater = greater
    def val(self):
        return self.val
    def inequality(self):
        return self.greater
class GPoint(Point):
    """
    Planar point
    """
    def __init__(self,pos,s_,d_=None):
        Point.__init__(self,pos)        
        if d_ is None:
            self.dir = s_
        else:
            self.strike = s_
            self.dip = d_
            self.dir = np.zeros(3)
            s_r = np.deg2rad(s_)
            d_r = np.deg2rad(np.abs(d_))
            self.dir[0] = m.sin(d_r)*m.cos(s_r)
            self.dir[1] = -m.sin(d_r)*m.sin(s_r)
            self.dir[2] = m.cos(d_r)
            self.dir  /= nla.norm(self.dir)
    def dir_(self):
        return self.dir
    def val(self,i):
        return self.dir[i]
class TPoint(Point):
    """
    Tangent point
    """
    def __init__(self,pos,s_,d_=None):
        Point.__init__(self,pos)        
        if d_ is None:
            self.dir = s_
        else:
            self.dir = np.zeros(3)
            s_r = np.deg2rad(s_)
            d_r = np.deg2rad(np.abs(d_))
            self.dir[0] = m.sin(d_r)*m.cos(s_r)
            self.dir[1] = -m.sin(d_r)*m.sin(s_r)
            self.dir[2] = m.cos(d_r)
            self.dir  /= nla.norm(self.dir)
    def dir_(self):
        return self.dir
    def val(self,i):
        return self.dir[i]
    def tx(self):
        return self.val(0)
    def ty(self):
        return self.val(1)
    def tz(self):
        return self.val(2)