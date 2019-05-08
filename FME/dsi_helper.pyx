import cython
import numpy as np
import numpy.linalg as la
from libc.math cimport sin, cos, M_PI, sqrt

cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
#cimport numpy.linalg as la
def get_element(long [:,:] elements, double [:,:] nodes, double [:] point):
    cdef int i
    cdef double[:, :] M = np.ones((4, 4))
    cdef double[:, :] points = np.empty((4, 3))#double [4][3] points
    cdef double[:] cp = np.ones(4)
    cdef double[:] c = np.empty(4)
    cdef int all_positive
    for i in range(3):
        cp[i+1] = point[i]
    for i in range(len(elements)):
        all_positive = 1
        for j in range(4):
            for k in range(3):
                M[j,k+1] = nodes[elements[i][j]][k]
        c = np.dot(cp,la.inv(M))
        for j in range(4):
            if c[j] < 0.: 
                all_positive = 0
                break
        if all_positive == 1:
            return elements[i]
    return None
def pointintetra(double[:,:] nodes, double [:] point):
    cdef int i
    cdef double[:, :] M = np.ones((4, 4))
    cdef double[:] cp = np.ones(4)
    cdef double[:] c = np.empty(4)
    cdef int all_positive
    for i in range(3):
        cp[i+1] = point[i]
    all_positive = 1
    for j in range(4):
        for k in range(3):
            M[j,k+1] = nodes[j][k]
    c = np.dot(cp,la.inv(M))
    for j in range(4):
        if c[j] < 0.: 
            all_positive = 0
            break
    if all_positive == 1:
        return True
    else:
        return False
def compute_regularisation_constraint(double [:,:] e1, double [:,:] e2, double [:] Xl, double [:] Xr, long [:] idl, long [:] idr):
    cdef int Nc, Na, i,position_to_write
    cdef long common_index, next_available_position
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    cdef double [:] c = np.zeros(Nc)
    cdef long [5] idc 
    for i in range(Nc):
        idc[i] = -1
    for itr_left in range(Na):
        idc[itr_left] = idl[itr_left]
        for i in range(3):
            c[itr_left] += Xl[i]*e1[i][itr_left]
    next_available_position = Na
    for itr_right in range(Na):
        common_index = -1
        for itr_left in range(Na):
            if idc[itr_left] == idr[itr_right]:
                common_index = itr_left

        position_to_write = 0
        if common_index != -1:
            position_to_write = common_index
        else:
            position_to_write = next_available_position
            next_available_position+=1
        idc[position_to_write] = idr[itr_right]
        for i in range(3):
            c[position_to_write] -= Xr[i]*e2[i][itr_right]
    return idc, c
def compute_gradient_ortho_constraint(double [:,:] element_gradient, double [:] direction ):
    cdef double [:] c = np.zeros(4)
    cdef int k, i
    for k in range(4):
        for i in range(3):
            c[k] += direction[i]*element_gradient[i][k]
    return c
def compute_cg_regularisation_constraint(double [:,:] e1,double [:,:] e2,long [:] idl, long [:] idr, double [:,:] nodes):
    cdef int Nc, Na, i,Ns, j
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    cdef double [:] c = np.zeros(Nc)
    cdef long [5] idc 
    cdef long [3] common 
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    for i in range(Nc):
        idc[i] = -1
    i = 0
    for itr_right in range(Na):
        for itr_left in range(Na):
            if idl[itr_left] == idr[itr_right]:
                common[i] = idl[itr_left]
                i+=1
    for j in range(i):
        for k in range(3):
            shared_pts[j][k] = nodes[common[j]][k]#common
    for i in range(3):
        v1[i] = shared_pts[0,i] - shared_pts[1,i]
        v2[i] = shared_pts[2,i]-shared_pts[1,i]
    norm = np.cross(v1,v2)
    #norm[0] = v1[1]*v2[2]-v1[2]*v2[1]
    #norm[1] = v1[2]*v2[0]-v1[0]*v2[2]
    #norm[2] = v1[0]*v2[1] - v1[1]*v2[0]
    for itr_left in range(Na):
        idc[itr_left] = idl[itr_left]
        for i in range(3):
            c[itr_left] += norm[i]*e1[i][itr_left]
    next_available_position = Na
    for itr_right in range(Na):
        common_index = -1
        for itr_left in range(Na):
            if idc[itr_left] == idr[itr_right]:
                common_index = itr_left

        position_to_write = 0
        if common_index != -1:
            position_to_write = common_index
        else:
            position_to_write = 4#next_available_position
            next_available_position+=1
        idc[position_to_write] = idr[itr_right]
        for i in range(3):
            c[position_to_write] -= norm[i]*e2[i][itr_right]
    return idc, c
def rotation(double [:] axis, double angle, double [:] vector, double [:,:] rotation_mat):
    cdef double c, s, C, x, y, z, xs, ys, xC, yC, zC, xyC, yzC, zxC, n
    cdef int i, j
    cdef double [:] rv = np.zeros(3)
    c = cos(angle*M_PI/180.)#np.cos(np.deg2rad(angle))
    s = sin(angle*M_PI/180.)#np.sin((np.deg2rad(angle)))
    C = 1.0 - c
    x = axis[0]
    y = axis[1]
    z = axis[2]
    xs = x*s
    ys = y*s
    zs = z*s
    xC = x*C
    yC = y*C
    zC = z*C
    xyC = x*yC
    yzC = y*zC
    zxC = z*xC
    rotation_mat[0][0] = x*xC+c
    rotation_mat[0][1] = xyC-zs
    rotation_mat[0][2] = zxC+ys
    
    rotation_mat[1][0] = xyC+zs
    rotation_mat[1][1] = y*yC+c
    rotation_mat[1][2] = yzC-xs

    rotation_mat[2][0] = zxC -ys
    rotation_mat[2][1] = yzC+xs
    rotation_mat[2][2] = z*zC+c
    
    for i in range(3):
        for j in range(3):
            rv[i]+=rotation_mat[i][j]*vector[i]
        n+=rv[i]*rv[i]
    n = sqrt(n)
    for i in range(3):
        rv[i]/=n
    return rv
def fourier_series(double [:] wl, double [:] C, double s):
    cdef double v
    cdef int i, w, wl_n, N
    N = 1
    wl_n = len(wl)
    for w in range(wl_n):
        for i in range(1,N+1):
             v += C[(2*i-1)+2*N*w]*cos(2*M_PI/wl[w] * i * (s)) + \
                C[(2*i)+2*N*w]*sin(2*M_PI/wl[w] * i * (s))
    return v
def cg(double [:,:,:] EG, long [:,:] neighbours, long [:,:] elements,double [:,:] nodes):
    cdef int Nc, Na, i,Ns, j, ne, ncons, e, n, neigh
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    ne = len(neighbours)
    ncons = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef double [:,:] c = np.zeros((len(neighbours)*4,Nc))
    cdef long [:,:] idc = np.zeros((ne*4,5),dtype=np.int64)
    cdef long [3] common 
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    cdef double [:,:] e1
    cdef double [:,:] e2
    cdef long [:] idl  = np.zeros(4,dtype=np.int64) 
    cdef long [:] idr = np.zeros(4,dtype=np.int64) 
    for e in range(ne):
        idl = elements[e,:]
        e1 = EG[e,:,:]
        flag[e] = 1
        for n in range(4):
            
            neigh = neighbours[e,n]
            idr = elements[neigh,:]
            if flag[neigh]== 1:
                continue
            if neigh == -1:
                continue
            e2 = EG[neigh,:,:]


            
            for i in range(Nc):
                idc[ncons,i] = -1

            i = 0
            for itr_right in range(Na):
                for itr_left in range(Na):
                    if idl[itr_left] == idr[itr_right]:
                        common[i] = idl[itr_left]
                        i+=1
            for j in range(3):
                for k in range(3):
                    shared_pts[j][k] = nodes[common[j]][k]#common
            for i in range(3):
                v1[i] = shared_pts[0,i] - shared_pts[1,i]
                v2[i] = shared_pts[2,i]-shared_pts[1,i]
            norm[0] = v2[2]*v1[1] - v1[2]*v2[1]
            norm[1] = v1[2]*v2[0] - v1[0]*v2[2]
            norm[2] = v1[0]*v2[1] - v1[1]*v2[0]#= np.cross(v1,v2)
            #norm[0] = v1[1]*v2[2]-v1[2]*v2[1]
            #norm[1] = v1[2]*v2[0]-v1[0]*v2[2]
            #norm[2] = v1[0]*v2[1] - v1[1]*v2[0]
            for itr_left in range(Na):
                idc[ncons,itr_left] = idl[itr_left]
                for i in range(3):
                    c[ncons,itr_left] += norm[i]*e1[i][itr_left]
            next_available_position = Na
            for itr_right in range(Na):
                common_index = -1
                for itr_left in range(Na):
                    if idc[ncons,itr_left] == idr[itr_right]:
                        common_index = itr_left

                position_to_write = 0
                if common_index != -1:
                    position_to_write = common_index
                else:
                    position_to_write = 4#next_available_position
                    next_available_position+=1
                idc[ncons,position_to_write] = idr[itr_right]
                for i in range(3):
                    c[ncons,position_to_write] -= norm[i]*e2[i][itr_right]
            ncons+=1
    return idc, c, ncons
def fold_cg(double [:,:,:] EG, double [:,:] X, long [:,:] neighbours, long [:,:] elements,double [:,:] nodes):
    cdef int Nc, Na, i,Ns, j, ne, ncons, e, n, neigh
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    ne = len(neighbours)
    ncons = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef double [:,:] c = np.zeros((len(neighbours)*4,Nc))
    cdef long [:,:] idc = np.zeros((ne*4,5),dtype=np.int64)
    cdef long [3] common 
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    cdef double [:,:] e1
    cdef double [:,:] e2
    cdef double [:] Xl
    cdef double [:] Xr

    cdef long [:] idl  = np.zeros(4,dtype=np.int64) 
    cdef long [:] idr = np.zeros(4,dtype=np.int64) 
    for e in range(ne):
        idl = elements[e,:]
        e1 = EG[e,:,:]
        flag[e] = 1
        Xl = X[e,:]
        for n in range(4):
            neigh = neighbours[e,n]
            idr = elements[neigh,:]
            if flag[neigh]== 1:
                continue
            if neigh == -1:
                continue
            e2 = EG[neigh,:,:]
            Xr = X[neigh,:]


            
            for i in range(Nc):
                idc[ncons,i] = -1

            i = 0
            for itr_left in range(Na):
                idc[ncons,itr_left] = idl[itr_left]
                for i in range(3):
                    c[ncons,itr_left] += Xl[i]*e1[i][itr_left]
            next_available_position = Na
            for itr_right in range(Na):
                common_index = -1
                for itr_left in range(Na):
                    if idc[ncons,itr_left] == idr[itr_right]:
                        common_index = itr_left
                position_to_write = 0
                if common_index != -1:
                    position_to_write = common_index
                else:
                    position_to_write = next_available_position
                    next_available_position+=1
                idc[ncons,position_to_write] = idr[itr_right]
                for i in range(3):
                    c[ncons,position_to_write] -= Xr[i]*e2[i][itr_right]
            ncons+=1
    return idc, c, ncons
def cg_cstr_to_coo_sq(double [:,:] c, long [:,:] idc, long ncons):
    cdef double [:] A = np.zeros(ncons*5*5)
    cdef int nc
    cdef long [:] row = np.zeros(ncons*5*5).astype(np.int64)
    cdef long [:] col = np.zeros(ncons*5*5).astype(np.int64)
    cdef long counter
    counter = 0
    nc = c.shape[1]
    for i in range(ncons):
        for j in range(nc):
            for l in range(nc):
                A[counter] = c[i,j]*c[i,l]
                row[counter] = idc[i,l]
                col[counter] = idc[i,j]
                counter+=1
    return A, row, col
def cg_cstr_to_coo_sq_B(double [:,:] c, double [:] Bc, long [:,:] idc, long ncons, long nx):
    cdef double [:] A = np.zeros(ncons*5*5)
    cdef double [:] B = np.zeros(nx)
    cdef int nc
    cdef long [:] row = np.zeros(ncons*5*5).astype(np.int64)
    cdef long [:] col = np.zeros(ncons*5*5).astype(np.int64)
    cdef long counter
    counter = 0
    nc = c.shape[1]
    for i in range(ncons):
        for j in range(nc):
            B[idc[i,j]]+=c[i,j]*Bc[i]
            for l in range(nc):
                A[counter] = c[i,j]*c[i,l]
                row[counter] = idc[i,l]
                col[counter] = idc[i,j]
                counter+=1
    return A, B, row, col