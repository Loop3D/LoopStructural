import numpy as np
import math as m
from numpy import linalg as nla


class Point():
    """
    Base object for a point contains the geospatial location
    """
    def __init__(self, pos):
        self.type = 'Point'
        self.pos = np.array(pos)
        self.orig = np.array(pos)
    def transform(self,t):
        t /= nla.norm(t)
        t[0] = - t[1]
        t[1] = t[0]
        t[2] = t[2]
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
        if point.pos.ndim == 2:
            raise BaseException
        d = 0.0
        d = dist.dist_1d(point.pos,self.pos,d)+eps
        return d,point.pos-self.pos
    def dim(self):
        return len(self.pos)

class IPoint(Point):
    """
    Interface point
    """
    def __init__(self,pos,val):
        Point.__init__(self,pos)
        self.val = val
        self.type = 'IPoint'
    def val(self):
        return self.val

class IePoint(Point):
    """
    Inequality point
    """
    def __init__(self,pos,val,greater=True):
        Point.__init__(self,pos)
        self.type = 'IePoint'
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
    def __init__(self,pos,vec):
        Point.__init__(self,pos)
        self.type = 'GPoint'
        self.vec = vec
    @classmethod
    def from_plunge_plunge_dir(cls, pos, plunge,plunge_dir, polarity=1):
        plunge = np.deg2rad(plunge)
        plunge_dir = np.deg2rad(plunge_dir+90)
        vec = np.zeros(3)
        vec[0] = m.sin(plunge) * m.cos(plunge_dir)
        vec[1] = -m.sin(plunge) * m.sin(plunge_dir)
        vec[2] = m.cos(plunge)
        vec /= nla.norm(vec)
        vec*=polarity
        return cls(pos,vec)
    @classmethod
    def from_strike_and_dip(cls, pos, strike, dip, polarity=1):
        dir = np.zeros(3)
        strike = np.deg2rad(strike)
        dip = np.deg2rad(np.abs(dip))
        dir[0] = m.sin(dip) * m.cos(strike)
        dir[1] = -m.sin(dip) * m.sin(strike)
        dir[2] = m.cos(dip)
        dir /= nla.norm(dir)
        dir*=polarity
        return cls(pos, dir)
    @classmethod
    def from_dip_dip_dir(cls, pos, dip_dir, dip, polarity=1):
        dir = np.zeros(3)
        strike = np.deg2rad(dip_dir+90)
        dip = np.deg2rad(dip)
        dir[0] = m.sin(dip) * m.cos(strike)
        dir[1] = -m.sin(dip) * m.sin(strike)
        dir[2] = m.cos(dip)
        dir /= nla.norm(dir)
        dir*=polarity
        return cls(pos, dir)
    def dir_(self):
        return self.vec
    def val(self,i):
        return self.vec[i]
class TPoint(Point):
    """
    Tangent point
    """
    def __init__(self,pos,s_,d_=None):
        Point.__init__(self,pos)
        self.type = 'TPoint'
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
