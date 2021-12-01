import logging

import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class CubicFunction:
    """

    """
    def __init__(self):
        """
        Class to represent a cubic function
        """
        self.A = []  # np.zeros((4,4))
        self.B = []  # np.zeros((4))
        self.max_v = 999999
        self.min_v = -99999
        self.w = None

    def add_cstr(self, x, y):
        self.A.append([x ** 3, x ** 2, x, 1.])
        self.B.append(y)

    def add_grad(self, x, g):
        self.A.append([3 * x ** 2, 2 * x, 1., 0.])
        self.B.append(g)

    def add_max(self, max_v):
        self.max_v = max_v

    def add_min(self, min_v):
        self.min_v = min_v

    def __call__(self, v):
        if len(self.B) < 3:
            print("underdetermined")
            return
        if self.w is None:
            A = np.array(self.A)
            B = np.array(self.B)
            ATA = A.T @ A
            ATB = A.T @ B
            self.w = np.linalg.lstsq(ATA, ATB, rcond=None)[0]
        eva = self.w[0] * v ** 3 + self.w[1] * v ** 2 + self.w[2] * v + self.w[
            3]
        eva[v > self.max_v] = self.w[0] * self.max_v ** 3 + \
                              self.w[1] * self.max_v ** 2 + self.w[
                                  2] * self.max_v + self.w[3]
        eva[v < self.min_v] = self.w[0] * self.min_v ** 3 + \
                              self.w[1] * self.min_v ** 2 + self.w[
                                  2] * self.min_v + self.w[3]
        return eva


class Composite():
    """

    """
    def __init__(self, positive, negative):
        self.positive = positive
        self.negative = negative

    def __call__(self, v):
        v = np.array(v)
        r = np.zeros(v.shape)
        r[v > 0] = self.positive(v[v > 0])
        r[v < 0] = self.negative(v[v < 0])
        return r


class Ones:
    """

    """
    def __init__(self):
        pass

    def __call__(self, v):
        v = np.array(v)
        return np.ones(v.shape)


class Zeros:
    """

    """
    def __init__(self):
        pass

    def __call__(self, v):
        v = np.array(v)
        return np.zeros(v.shape)


class FaultDisplacement:
    """

    """
    def __init__(self, hw=None, fw=None, gx=None, gy=None, gz=None, **kwargs):
        self.gx = gx
        if hw is not None and fw is not None:
            self.__gx = Composite(hw, fw)
        self.gy = gy
        self.gz = gz
        self.gx_bounds = None
        self.gy_bounds = None
        self.gz_bounds = None

        if self.gx == None:
            print('Gx function none setting to ones')
            self.gx = Ones()
        if self.gy == None:
            print('Gy function none setting to ones')
            self.gy = Ones()
        if self.gz == None:
            print('Gz function none setting to ones')
            self.gz = Ones()
        if 'gxmax' in kwargs and 'gxmin' in kwargs:
            self.gx_bounds = (kwargs['gxmin'], kwargs['gxmax'])

        if 'gymax' in kwargs and 'gymin' in kwargs:
            self.gy_bounds = (kwargs['gymin'], kwargs['gymax'])

        if 'gzmax' in kwargs and 'gzmin' in kwargs:
            self.gz_bounds = (kwargs['gzmin'], kwargs['gzmax'])

    def __call__(self, gx, gy, gz):
        if self.gx_bounds is not None:
            mid =(self.gx_bounds[1]+self.gx_bounds[0])/2.
            gx = (gx -mid) / (self.gx_bounds[1]-self.gx_bounds[0])
        if self.gy_bounds is not None:
            mid =  (self.gy_bounds[1]+self.gy_bounds[0])/2.
            gy = (gy -mid) / (self.gy_bounds[1]-self.gy_bounds[0])
        if self.gz_bounds is not None:
            mid =(self.gz_bounds[1]+self.gz_bounds[0])/2.
            gz = (gz -mid) / (self.gz_bounds[1]-self.gz_bounds[0])
        return self.gx(gx)*self.gy(gy)*self.gz(gz)

    def evaluate(self, gy, gz):
        if self.gy_bounds is not None:
            mid =  (self.gy_bounds[1]+self.gy_bounds[0])/2.
            gy = (gy -mid) / (self.gy_bounds[1]-self.gy_bounds[0])
        if self.gz_bounds is not None:
            mid =(self.gz_bounds[1]+self.gz_bounds[0])/2.
            gz = (gz -mid) / (self.gz_bounds[1]-self.gz_bounds[0])
        return self.gy(gy)*self.gz(gz)


class BaseFault(object):
    """

    """
    hw = CubicFunction()
    hw.add_cstr(0, 1)
    hw.add_grad(0, 0)
    hw.add_cstr(1, 0)
    # hw.add_cstr(1,1)

    hw.add_grad(1, 0)
    hw.add_max(1)
    fw = CubicFunction()
    fw.add_cstr(0, -1)
    fw.add_grad(0, 0)
    fw.add_cstr(-1, 0)
    fw.add_grad(-1, 0)
    fw.add_min(-1)
    # gyf = CubicFunction()
    # gyf.add_cstr(-1, 0)
    # gyf.add_cstr(1, 0)
    # gyf.add_cstr(-0.2, 1)
    # gyf.add_cstr(0.2, 1)
    # gyf.add_grad(0, 0)
    # gyf.add_min(-1)
    # gyf.add_max(1)
    gyf = Ones()
    gzf = CubicFunction()
    gzf.add_cstr(-1, 0)
    gzf.add_cstr(1, 0)
    gzf.add_cstr(-0.2, 1)
    gzf.add_cstr(0.2, 1)
    gzf.add_grad(0, 0)
    gzf.add_min(-1)
    gzf.add_max(1)
    gxf = Composite(hw, fw)
    fault_displacement = FaultDisplacement(gx=gxf, gy=gyf, gz=gzf)