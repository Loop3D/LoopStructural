import numpy as np

class CubicFunction:
    def __init__(self):
        self.A = []#np.zeros((4,4))
        self.B = []#np.zeros((4))
        self.max_v = 999999
        self.min_v = -99999

    def add_cstr(self,x,y):
        self.A.append([x**3, x**2, x, 1.])
        self.B.append(y)

    def add_grad(self,x,g):
        self.A.append([3 * x ** 2, 2 * x, 1., 0.])
        self.B.append(g)

    def add_max(self, max_v):
        self.max_v = max_v

    def add_min(self, min_v):
        self.min_v = min_v

    def __call__(self,v):
        if len(self.B)<3:
            print("underdetermined")
            return
        A = np.array(self.A)
        B = np.array(self.B)
        ATA = A.T @ A
        ATB = A.T @ B
        w = np.linalg.lstsq(ATA,ATB)[0]
        eva = w[0]*v**3+w[1]*v**2+w[2]*v+w[3]
        eva[v>self.max_v] = w[0]*self.max_v**3+w[1]*self.max_v**2+w[2]*self.max_v+w[3]
        eva[v<self.min_v] =  w[0]*self.min_v**3+w[1]*self.min_v**2+w[2]*self.min_v+w[3]
        return eva

class Ones:
    def __init__(self):
        pass

    def __call__(self,v):
        v = np.array(v)
        return np.ones(v.shape)


class FaultDisplacement:
    def __init__(self,fw=None,hw=None,gy=None,gz=None):
        self.__fw = fw
        self.__hw = hw
        self.__gy = gy
        self.__gz = gz

    def fw(self,gx,gy,gz):
        return self.__fw(gx)*self.__gy(gy)*self.__gz(gz)

    def hw(self,gx,gy,gz):
        return self.__hw(gx)*self.__gy(gy)*self.__gz(gz)