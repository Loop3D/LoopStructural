import cython
import numpy as np
import numpy.linalg as la
cimport numpy as np
from math import *
@cython.boundscheck(False)
@cython.wraparound(False)
#cython: language_level=3
#cimport numpy.linalg as la
def cg(double [:,:,:] EG, long long [:,:] neighbours, long long [:,:] elements,double [:,:] nodes, long long [:] region):
    cdef int Nc, Na, i,Ns, j, ne, ncons, e, n, neigh
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    ne = len(neighbours)
    ncons = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef double [:,:] c = np.zeros((len(neighbours)*4,Nc))
    cdef double [:] areas = np.zeros((len(neighbours)*4))
    cdef long long [:,:] idc = np.zeros((ne*4,5),dtype=np.int64)
    cdef long long [3] common
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    cdef double [:,:] e1
    cdef double [:,:] e2
    cdef double area = 0
    cdef double length
    cdef long long [:] idl  = np.zeros(4,dtype=np.int64)
    cdef long long [:] idr = np.zeros(4,dtype=np.int64)
    for e in range(ne):
        idl = elements[e,:]
        e1 = EG[e,:,:]
        flag[e] = 1
        # if not in region then skip this tetra
        if region[idl[0]] == 0 or region[idl[1]] == 0 or region[idl[2]] == 0 or region[idl[3]] == 0:
            continue
        for n in range(4):
            neigh = neighbours[e,n]
            idr = elements[neigh,:]
            if neigh < 0:
                continue
            if flag[neigh]== 1:
                continue

            # if not in region then skip this tetra
            if region[idr[0]] == 0 or region[idr[1]] == 0 or region[idr[2]] == 0 or region[idr[3]] == 0:
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
            norm[2] = v1[0]*v2[1] - v1[1]*v2[0]

            length = np.linalg.norm(norm)
            # we want to weight the cg by the area of the shared face
            # area of triangle is half area of parallelogram
            # https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle
            area = 0.5*length#sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2])#np.linalg.norm(norm)
            for i in range(3):
                norm[i]/=length
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
            areas[ncons] = area
            ncons+=1
    return idc, c, ncons, areas
def constant_norm(double [:,:,:] EG, long long [:,:] neighbours, long long [:,:] elements,double [:,:] nodes, long long [:] region):
    cdef int Nc, Na, i,Ns, j, ne, ncons, e, n, neigh
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    ne = len(neighbours)
    ncons = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef double [:,:] c = np.zeros((len(neighbours)*4,Nc))
    cdef long long [:,:] idc = np.zeros((ne*4,5),dtype=np.int64)
    cdef long long [3] common
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    cdef double [:,:] e1
    cdef double [:,:] e2
    cdef double area = 0
    cdef long long [:] idl  = np.zeros(4,dtype=np.int64)
    cdef long long [:] idr = np.zeros(4,dtype=np.int64)
    for e in range(ne):
        idl = elements[e,:]
        e1 = EG[e,:,:]
        flag[e] = 1
        # if not in region then skip this tetra
        if region[idl[0]] == 0 or region[idl[1]] == 0 or region[idl[2]] == 0 or region[idl[3]] == 0:
            continue
        for n in range(4):
            neigh = neighbours[e,n]
            idr = elements[neigh,:]
            if neigh < 0:
                continue
            if flag[neigh]== 1:
                continue

            # if not in region then skip this tetra
            if region[idr[0]] == 0 or region[idr[1]] == 0 or region[idr[2]] == 0 or region[idr[3]] == 0:
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
            norm[2] = v1[0]*v2[1] - v1[1]*v2[0]

            # we want to weight the cg by the area of the shared face
            # area of triangle is half area of parallelogram
            # https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle
            area = 0.5*np.linalg.norm(norm)
            for itr_left in range(Na):
                idc[ncons,itr_left] = idl[itr_left]
                c[ncons,itr_left] = np.sqrt(e1[0][itr_left]*e1[0][itr_left]+e1[1][itr_left]*e1[1][itr_left]+e1[2][itr_left]*e1[2][itr_left])
              
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
                c[ncons,position_to_write] -= np.sqrt(e2[0][itr_right]*e2[0][itr_right]+e2[1][itr_right]*e2[1][itr_right]+e2[2][itr_right]*e2[2][itr_right]) 
            ncons+=1
    return idc, c, ncons
def fold_cg(double [:,:,:] EG, double [:,:] X, long long [:,:] neighbours, long long [:,:] elements,double [:,:] nodes):
    cdef int Nc, Na, i,Ns, j, ne, ncons, e, n, neigh
    Nc = 5 #numer of constraints shared nodes + independent
    Na = 4 #number of nodes
    Ns = Na -1
    ne = len(neighbours)
    ncons = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef double [:,:] c = np.zeros((len(neighbours)*4,Nc))
    cdef long long [:,:] idc = np.zeros((ne*4,5),dtype=np.int64)
    cdef long long [3] common
    cdef double [:] norm = np.zeros((3))
    cdef double [:,:] shared_pts = np.zeros((3,3))
    cdef double [:] v1 = np.zeros(3)
    cdef double [:] v2 = np.zeros(3)
    cdef double [:,:] e1
    cdef double [:,:] e2
    cdef double [:] Xl
    cdef double [:] Xr

    cdef long long [:] idl  = np.zeros(4,dtype=np.int64)
    cdef long long [:] idr = np.zeros(4,dtype=np.int64)
    for e in range(ne):
        idl = elements[e,:]
        e1 = EG[e,:,:]
        flag[e] = 1
        Xl = X[e,:]
        for n in range(4):
            neigh = neighbours[e,n]
            idr = elements[neigh,:]
            if neigh == -1:
                continue
            if flag[neigh]== 1:
                continue
            e2 = EG[neigh,:,:]
            Xr = X[neigh,:]



            for i in range(Nc):
                idc[ncons,i] = -1
            i=0
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
            norm[2] = v1[0]*v2[1] - v1[1]*v2[0]
            area = 0.5*sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2])#np.linalg.norm(norm)

            i = 0
            for itr_left in range(Na):
                idc[ncons,itr_left] = idl[itr_left]
                for i in range(3):
                    c[ncons,itr_left] += Xl[i]*e1[i][itr_left]*area
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
                    c[ncons,position_to_write] -= Xr[i]*e2[i][itr_right]*area
            ncons+=1
    return idc, c, ncons

def tetra_neighbours(long long [:,:] elements, long long [:,:] neighbours):
    cdef long long ie, ne, nn, n, i, j
    for ie in range(len(elements)):
        nn = 0 ## counter for number of neighbours
        for ne in range(len(elements)):
            n = 0 # counter for how many shared nodes
            if ne == ie:
                continue
            for i in range(4):
                for j in range(4):
                    if elements[ie,i] == elements[ne,j]:
                        n+=1
            if n == 3:
                neighbours[ie,nn] = ne
                nn+=1

def calculate_pairs(neighbours,elements):
    ne = len(neighbours)
    cdef long long [:,:] faces = np.zeros((ne*4,3),dtype=np.int64)
    cdef long long [:,:] pairs = np.zeros((ne*4,2),dtype=np.int64)    
    face_n = 0
    cdef int [:] flag = np.zeros(ne,dtype=np.int32)
    cdef int e = 0
    cdef int n =0
    for e in range(ne):
        flag[e] = 1
        idl = elements[e,:]
        for n in range(4):
            neigh = neighbours[e][n]
            idr = elements[neigh,:]
            if neigh < 0:
                continue
            if flag[neigh]== 1:
                continue
            pairs[face_n][0] = e
            pairs[face_n][1] = n

            face_n+=1
    return pairs