cimport numpy as np
import cython

from libc.math cimport sqrt,pow
#DTYPE = np.float 
#ctypedef np.float_t DTYPE_t 
@cython.boundscheck(False)
@cython.wraparound(False)
def dist_1d(double [:] points, double [:] p, double d):
    cdef int dim = len(points)#.shape[1]
    #cdef double d[1][:]  
    cdef int i
    cdef double t = 0.0

    for i in range(dim):
        t = t+pow((points[i]-p[i]),2)
    if t == 0.0:
        return 0.0
    d = sqrt(t)
    return d
def dist_2d(double [:,:] points, double [:] p, double [:] d, double [:,:] delta):
    cdef int npoints = len(points)#.shape[0]
    cdef int dim = len(points[1])#.shape[1]
    #cdef double d[1][:]  
    cdef int i, j
    cdef double t = 0.0
    for j in range(npoints):
        t = 0.0
        for i in range(dim):
            t = t+pow((points[j,i]-p[i]),2)
            delta[j,i] = points[j,i]-p[i]
        if t == 0.0:
            d[j] = 0.0
        else:
            d[j] = sqrt(t)
    return True #trued, delta

