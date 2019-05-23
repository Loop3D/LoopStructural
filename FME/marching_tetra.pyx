#cython: boundscheck=False, wraparound=False, nonecheck=False, linetrace=True
import numpy as np
cdef node_exists(nodes,node,n):
    """
    iterate through a container of nodes to see if the node exists if it exists return the index
    """
    cdef double csum
    if n == 0:
        return -1
    for i in range(n):
        csum = 0
        for j in range(3):
            csum += (nodes[i,j]-node[j])**2
        if csum == 0:
            return i
    return -1
cdef vertex_interp(double isolevel,double [:] v1, double [:] v2,double [:] point):
    """
    Find the coordinates of vertices on an edge where the isovalue intersects
    isovalue is a float of
    """
    cdef int i
    cdef double[:] p1, p2
    cdef double vp1, vp2, mu
    p1 = v1[:3]
    p2 = v2[:3]
    vp1 = v1[3]
    vp2 = v2[3]
    #print(vp1,vp2,isolevel)
    if abs(isolevel-vp1) < 0.00001 : 
        for i in range(3):
            point[i] = p1[i]
        return  
    if abs(isolevel-vp2) < 0.00001 : 
        for i in range(3):
            point[i] = p2[i]
        return
    ##if the value of both points is really similar return p1
    if abs(vp1-vp2) < 0.00001 : 
        for i in range(3):
            point[i] = p1[i]
        return  
    ## interpolate on a line
    mu = (isolevel - vp1) / (vp2 - vp1)  
    for i in range(3):
        point[i]= p1[i] + mu * (p2[i] - p1[i])
                
cdef t0E01(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):
    vertex_interp(iso,v0,v1,triangles[ntri[0],0,:])        
    vertex_interp(iso,v0,v2,triangles[ntri[0],1,:])
    vertex_interp(iso,v0,v3,triangles[ntri[0],2,:])
    ntri[0]+=1
    return 
  
cdef t0D02(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):  
    vertex_interp(iso,v1,v0,triangles[ntri[0],0,:])#g.p[v1],g.p[v0],g.val[v1],g.val[v0])
    vertex_interp(iso,v1,v3,triangles[ntri[0],1,:])#g.p[v1],g.p[v3],g.val[v1],g.val[v3])
    vertex_interp(iso,v1,v2,triangles[ntri[0],2,:])#g.p[v1],g.p[v2],g.val[v1],g.val[v2])
    ntri[0]+=1
    return 
  
cdef t0C03(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):  
    vertex_interp(iso,v0,v3,triangles[ntri[0],0,:])#g.p[v0],g.p[v3],g.val[v0],g.val[v3])
    vertex_interp(iso,v0,v2,triangles[ntri[0],1,:])#g.p[v0],g.p[v2],g.val[v0],g.val[v2])
    vertex_interp(iso,v1,v3,triangles[ntri[0],2,:])#g.p[v1],g.p[v3],g.val[v1],g.val[v3])
    ntri[0]+=1
    vertex_interp(iso,v1,v2,triangles[ntri[0],1,:])#g.p[v1],g.p[v2],g.val[v1],g.val[v2])
    for i in range(3):
        triangles[ntri[0],0,i] = triangles[ntri[0]-1,2,i]#g.p[v1],g.p[v3],g.val[v1],g.val[v3])
        triangles[ntri[0],2,i] = triangles[ntri[0]-1,1,i]#g.p[v0],g.p[v2],g.val[v0],g.val[v2]) 
    ntri[0]+=1
    return
  
cdef t0B04(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):
    vertex_interp(iso,v2,v0,triangles[ntri[0],0,:])#g.p[v2],g.p[v0],g.val[v2],g.val[v0])
    vertex_interp(iso,v2,v1,triangles[ntri[0],1,:])#g.p[v2],g.p[v1],g.val[v2],g.val[v1])
    vertex_interp(iso,v2,v3,triangles[ntri[0],2,:])#g.p[v2],g.p[v3],g.val[v2],g.val[v3])
    ntri[0]+=1
    return
cdef t0A05(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3): 
    vertex_interp(iso,v0,v1,    triangles[ntri[0],0,:])#g.p[v0],g.p[v1],g.val[v0],g.val[v1])
    vertex_interp(iso,v2,v3,triangles[ntri[0],1,:])#g.p[v2],g.p[v3],g.val[v2],g.val[v3])
    vertex_interp(iso,v0,v3,triangles[ntri[0],2,:])#g.p[v0],g.p[v3],g.val[v0],g.val[v3])
    ntri[0]+=1
    vertex_interp(iso,v1,v2,triangles[ntri[0],1,:])#g.p[v1],g.p[v2],g.val[v1],g.val[v2])
    for i in range(3):
        triangles[ntri[0],0,i] = triangles[ntri[0]-1,0,i]
        triangles[ntri[0],2,i] = triangles[ntri[0]-1,1,i]
    ntri[0]+=1
    return
  
cdef t0906(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):
    vertex_interp(iso,v0,v1,triangles[ntri[0],0,:])#g.p[v0],g.p[v1],g.val[v0],g.val[v1])
    vertex_interp(iso,v1,v3,triangles[ntri[0],1,:])#g.p[v1],g.p[v3],g.val[v1],g.val[v3])
    vertex_interp(iso,v2,v3,triangles[ntri[0],2,:])#g.p[v2],g.p[v3],g.val[v2],g.val[v3])
    ntri[0]+=1
    for i in range(3):
        triangles[ntri[0],0,i] = triangles[ntri[0]-1,0,i]
        triangles[ntri[0],2,i] = triangles[ntri[0]-1,2,i]
    vertex_interp(iso,v0,v2,triangles[ntri[0],1,:])#g.p[v0],g.p[v2],g.val[v0],g.val[v2])
    ntri[0]+=1
    return
cdef t0708(double [:,:,:] triangles, long[:] ntri,double iso, double [:] v0, double [:] v1, double [:] v2, double [:] v3):
    vertex_interp(iso,v3,v0,triangles[ntri[0],0,:])#g.p[v3],g.p[v0],g.val[v3],g.val[v0])
    vertex_interp(iso,v3,v2,triangles[ntri[0],1,:])#g.p[v3],g.p[v2],g.val[v3],g.val[v2])
    vertex_interp(iso,v3,v1,triangles[ntri[0],2,:])#g.p[v3],g.p[v1],g.val[v3],g.val[v1])
    ntri[0]+=1
    return 
def marching_tetra(double isovalue,long [:,:] elements,double [:,:] nodes, propertyvalue):
    nodes = np.hstack([nodes,propertyvalue[:,None]]) 
    #find which nodes are > isovalue
    property_bool = propertyvalue > isovalue
    #find what case each tetra is by slicing property bool by elements and multiplying by the case array
    #the sum of these rows = the case
    tetra_type_index1 = np.sum(property_bool[elements]*np.array([1,2,4,8]),axis=1)
    #0 and 15 are tcases where isovalue doesn't intersect surface so just skip these 
    tetra_type_index = tetra_type_index1[np.logical_and(tetra_type_index1>0,tetra_type_index1<15)]
    tetras_index = np.array(range(0,len(elements)))[np.logical_and(tetra_type_index1>0,tetra_type_index1<15)]
    ##create container for the triangle nodes and indexes and a counter
    cdef int ne, nn, i
    #use np array to make use of buffer
    cdef long[:] ntri = np.zeros(1).astype(int)
    cdef double[:] v0,v1,v2,v3
    ne =0
    nn =0
    ntri[0]=0
    cdef double [:,:,:] triangles = np.zeros((len(nodes),3,3))
    for i in range(len(tetra_type_index)):
        #get only the nodes of the tetras that containt the surface
        #vertex is x,y,z,propertyval
        v0 = nodes[elements[tetras_index[i]][0]]
        v1 = nodes[elements[tetras_index[i]][1]]
        v2 = nodes[elements[tetras_index[i]][2]]
        v3 = nodes[elements[tetras_index[i]][3]]     
        if tetra_type_index[i]== 7:
            t0708(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]== 8:
            t0708(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==9:
            t0906(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==6:
            t0906(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]== 10:
            t0A05(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==5:
            t0A05(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==11:
            t0B04(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]== 4:
            t0B04(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==12:
            t0C03(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==3:
            t0C03(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==13:
            t0D02(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]== 2:
            t0D02(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==14:
            t0E01(triangles,ntri,isovalue,v0,v1,v2,v3)
        if tetra_type_index[i]==1:
            t0E01(triangles,ntri,isovalue,v0,v1,v2,v3)
    return triangles, ntri
    #cdef double [:,:] triangle_nodes = np.zeros(((ntri[0]),3))
    #cdef long [:,:] triangle_indexes = np.zeros(((ntri[0]),3)).astype(int)
    #nn = 0 
    #for i in range(ntri[0]):
    #    for j in range(3):
    #        idx = node_exists(triangle_nodes,triangles[i,j,:],nn)
    #        if idx == -1:
    #            triangle_nodes[nn,:]=triangles[i,j,:]
    #            triangle_indexes[i,j] = nn
    #            nn+=1
    #        else:
    #            triangle_indexes[i,j] = idx
    #return triangle_nodes, triangle_indexes, nn
