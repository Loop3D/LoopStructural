"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np
from LoopStructural.interpolators.cython.dsi_helper import cg, constant_norm, fold_cg

logger = logging.getLogger(__name__)

class TriMesh:
    """

    """
    def __init__(self, elements, vertices,neighbours):
        self.elements = elements
        self.vertices = vertices
        if self.elements.shape[1] == 3:
            self.order =  1
        elif self.elements.shape[1] == 6:
            self.order = 2
        self.nelements = self.elements.shape[0]
        self.nx = self.vertices.shape[0]
        self.n_nodes  = self.nx
        self.neighbours = neighbours
        
        self.hN = np.array([[4,4,0,0,0,-8],[4,0,4,0,-8,0]])
        self.Nst = np.array([4,0,0,4,-4,-4])
        self.properties = {}
        # build an array of edges and edge relationships
        self.edges = np.zeros((self.nelements*3,2),dtype=int)
        self.edge_relationships = np.zeros((self.nelements*3,2),dtype=int)
        edge_index=0
        flag = np.zeros(self.elements.shape[0],dtype=bool)
        for i, t in enumerate(self.elements):
            flag[i] = True
            for n in self.neighbours[i]:
                if n < 0:
                    continue
                if flag[n]:
                    continue
                edge_node_index=0
                self.edge_relationships[edge_index,0] = i
                self.edge_relationships[edge_index,1] = n
                for v in t:
                    if v in self.elements[n,:3]:
                        self.edges[edge_index,edge_node_index] = v
                        edge_node_index+=1

                edge_index+=1
        self.edges = self.edges[:edge_index,:]
        self.edge_relationships = self.edge_relationships[:edge_index,:]
    @property
    def nodes(self):
        """
        Gets the nodes of the mesh as a property rather than using a function, accessible as a property! Python magic!

        Returns
        -------
        nodes : np.array((N,3))
            Fortran ordered
        """
        return self.vertices
    
    def barycentre(self, elements = None):
        """
        Return the barycentres of all tetrahedrons or of specified tetras using
        global index

        Parameters
        ----------
        elements - numpy array
            global index

        Returns
        -------

        """
        if elements is None:
            element_idx = np.arange(0,self.nelements)
        elements = self.get_elements()[element_idx]
        barycentre = np.sum(self.nodes[elements][:, :, :],
                                 axis=1) / 4.
        return barycentre

    def element_area(self, elements):
        tri_points = self.nodes[self.elements[elements,:],:]
        M_t = np.ones((tri_points.shape[0],3,3))
        M_t[:,:,1:] = tri_points[:,:3,:]    
        area = np.abs(np.linalg.det(M_t))*0.5
        return area
    def evaluate_value(self, pos, prop):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - numpy array
            property values at nodes

        Returns
        -------

        """
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        c, tri = self.evaluate_shape(pos[:,:2])
        inside = tri > 0
        # vertices, c, elements, inside = self.get_elements_for_location(pos)
        values[inside] = np.sum(c[inside,:]*self.properties[prop][self.elements[tri[inside],:]],axis=1)
        return values

    def evaluate_gradient(self, pos, prop):
        """
        Evaluate the gradient of an interpolant at the locations

        Parameters
        ----------
        pos - numpy array
            locations
        prop - string
            property to evaluate


        Returns
        -------

        """
        values = np.zeros(pos.shape)
        values[:] = np.nan
        element_gradients, tri = self.evaluate_shape_derivatives(pos[:,:2])
        inside = tri > 0
        # ?vertices, element_gradients, elements, inside = self.get_element_gradient_for_location(pos[:,:2])
        # vertex_vals = self.properties[prop][elements]
        #grads = np.zeros(tetras.shape)
        # v = (element_gradients[inside,:,:]*self.support.properties['defaultproperty'][self.support.elements[tri[inside],:,None]]).sum(1)
        values[inside,:] = (element_gradients[inside,:,:]*self.properties[prop][self.elements[tri[inside],:,None]]).sum(1)
        length = np.sum(values[inside,:],axis=1)
        # values[inside,:] /= length[:,None]
        return values


    def get_element_for_location(self, pos):
        """
        Determine the elements from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        M = np.ones((self.elements.shape[0],3,3))
        M[:,:,1:] = self.vertices[self.elements,:][:,:3,:]
        points_ = np.ones((pos.shape[0],3))
        points_[:,1:] = pos
        minv = np.linalg.inv(M)
        c = np.einsum('kij,li->lkj',minv,points_)
        isinside = np.all(c > 0,axis=2)
        ix,iy = np.where(isinside==True)
        element_idx = np.zeros(pos.shape[0],dtype=int)
        element_idx[:] = -1
        element_idx[ix] = iy 
        c_return = np.zeros((pos.shape[0],3))
        c_return[:] = np.nan
        c_return[ix,:] = c[isinside,:]
        return  c_return, element_idx

    def evaluate_mixed_derivative(self,indexes): 
        """
        evaluate partial of N with respect to st (to set u_xy=0)
        """

        vertices = self.nodes[self.elements[indexes],:]
        jac  = np.array([[(vertices[:,1,0]-vertices[:,0,0]),(vertices[:,1,1]-vertices[:,0,1])],
            [vertices[:,2,0]-vertices[:,0,0],vertices[:,2,1]-vertices[:,0,1]]]).T
        Nst_coeff = jac[:,0,0]*jac[:,1,1]+jac[:,0,1]*jac[:,1,0]

        
        Nst = self.Nst[None,:]*Nst_coeff[:,None]
        return Nst + self.hN[None,0,:] * (jac[:,0,0]*jac[:,1,0])[:,None] + self.hN[None,1,:] * (jac[:,1,0]*jac[:,1,1])[:,None]
  
    def evaluate_shape_d2(self,indexes): 
        """evaluate second derivatives of shape functions in s and t

        Parameters
        ----------
        M : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """

        vertices = self.nodes[self.elements[indexes],:]
        
        jac  = np.array([[(vertices[:,1,0]-vertices[:,0,0]),(vertices[:,1,1]-vertices[:,0,1])],
            [vertices[:,2,0]-vertices[:,0,0],vertices[:,2,1]-vertices[:,0,1]]]).T
        jac = np.linalg.inv(jac)
        jac = jac*jac


        
        d2_prod = np.einsum('lij,ik->lik',jac,self.hN)
        d2Const = d2_prod[:,0,:] + d2_prod[:,1,:]
        xxConst = d2_prod[:,0,:]
        yyConst = d2_prod[:,1,:]

        return xxConst,yyConst

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        compute dN/ds (1st row), dN/dt(2nd row)
        """
        locations = np.array(locations)
        if elements is None:
            c, tri =  self.get_element_for_location(locations)
        else:
            tri = elements
            M = np.ones((elements.shape[0],3,3))
            M[:,:,1:] = self.vertices[self.elements[elements],:][:,:3,:]
            points_ = np.ones((locations.shape[0],3))
            points_[:,1:] = locations
            minv = np.linalg.inv(M)
            c = np.einsum('lij,li->lj',minv,points_)

        vertices = self.nodes[self.elements[tri][:,:3]]
        jac = np.zeros((tri.shape[0],2,2))
        jac[:,0,0] = vertices[:,1,0]-vertices[:,0,0]
        jac[:,0,1] = vertices[:,1,1]-vertices[:,0,1]
        jac[:,1,0] = vertices[:,2,0]-vertices[:,0,0]
        jac[:,1,1] = vertices[:,2,1]-vertices[:,0,1]
        N = np.zeros((tri.shape[0],6))
        
        # dN containts the derivatives of the shape functions
        dN = np.zeros((tri.shape[0],2,6))
        dN[:,0,0] = 4*c[:,1]+4*c[:,2]-3#diff(N1,s).evalf(subs=vmap)
        dN[:,0,1] = 4*c[:,1]-1#diff(N2,s).evalf(subs=vmap)
        dN[:,0,2] = 0#diff(N3,s).evalf(subs=vmap)
        dN[:,0,3] = 4*c[:,2] #diff(N4,s).evalf(subs=vmap)
        dN[:,0,4] = -4*c[:,2]#diff(N5,s).evalf(subs=vmap)
        dN[:,0,5] = -8*c[:,1]-4*c[:,2]+4#diff(N6,s).evalf(subs=vmap)
        
        dN[:,1,0] = 4*c[:,1]+4*c[:,2]-3#diff(N1,t).evalf(subs=vmap)
        dN[:,1,1] = 0#diff(N2,t).evalf(subs=vmap)
        dN[:,1,2] = 4*c[:,2]-1#diff(N3,t).evalf(subs=vmap)
        dN[:,1,3] = 4*c[:,1]#diff(N4,t).evalf(subs=vmap)
        dN[:,1,4] = -4*c[:,1]-8*c[:,2]+4#diff(N5,t).evalf(subs=vmap)
        dN[:,1,5] = -4*c[:,1]#diff(N6,t).evalf(subs=vmap)
        
        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
#         d_n = np.einsum('ijk,ijl->ilk',np.linalg.inv(jac),dN)
        d_n = np.linalg.inv(jac)
#         d_n = d_n.swapaxes(1,2)
        d_n = d_n @ dN
        d_n = d_n.swapaxes(2,1)
        # d_n = np.dot(np.linalg.inv(jac),dN)
        return d_n, tri

    def evaluate_shape_derivativeso(self,x,y,M):
        """
        compute dN/ds (1st row), dN/dt(2nd row)
        """
        c = np.dot(np.array([1,x,y]),np.linalg.inv(M))
 
        vertices = np.zeros((3,2))
        vertices[0,:] = [M[0,1],M[0,2]]
        vertices[1,:] = [M[1,1],M[1,2]]
        vertices[2,:] = [M[2,1],M[2,2]]
        jac  = np.array([[(vertices[1,0]-vertices[0,0]),(vertices[1,1]-vertices[0,1])],
            [vertices[2,0]-vertices[0,0],vertices[2,1]-vertices[0,1]]])
        N = np.zeros(6)
        
        # dN containts the derivatives of the shape functions
        dN = np.zeros((2,6))
        dN[0,0] = 4*c[1]+4*c[2]-3#diff(N1,s).evalf(subs=vmap)
        dN[0,1] = 4*c[1]-1#diff(N2,s).evalf(subs=vmap)
        dN[0,2] = 0#diff(N3,s).evalf(subs=vmap)
        dN[0,3] = 4*c[2] #diff(N4,s).evalf(subs=vmap)
        dN[0,4] = -4*c[2]#diff(N5,s).evalf(subs=vmap)
        dN[0,5] = -8*c[1]-4*c[2]+4#diff(N6,s).evalf(subs=vmap)
        
        dN[1,0] = 4*c[1]+4*c[2]-3#diff(N1,t).evalf(subs=vmap)
        dN[1,1] = 0#diff(N2,t).evalf(subs=vmap)
        dN[1,2] = 4*c[2]-1#diff(N3,t).evalf(subs=vmap)
        dN[1,3] = 4*c[1]#diff(N4,t).evalf(subs=vmap)
        dN[1,4] = -4*c[1]-8*c[2]+4#diff(N5,t).evalf(subs=vmap)
        dN[1,5] = -4*c[1]#diff(N6,t).evalf(subs=vmap)
        
    
        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
        d_n = np.dot(np.linalg.inv(jac),dN)
    
        return d_n 

    def evaluate_shape(self,locations):
        locations = np.array(locations)
        c, tri = self.get_element_for_location(locations)
        #c = np.dot(np.array([1,x,y]),np.linalg.inv(M)) # convert to barycentric coordinates
        # order of bary coord is (1-s-t,s,t)
        N = np.zeros((c.shape[0],6)) #evaluate shape functions at barycentric coordinates
        N[:,0] = c[:,0]*(2*c[:,0]-1) #(1-s-t)(1-2s-2t)
        N[:,1] = c[:,1]*(2*c[:,1]-1) #s(2s-1)
        N[:,2] = c[:,2]*(2*c[:,2]-1) #t(2t-1)
        N[:,3] = 4*c[:,1]*c[:,2] #4st 
        N[:,4] = 4*c[:,2]*c[:,0] #4t(1-s-t)
        N[:,5] = 4*c[:,1]*c[:,0] #4s(1-s-t)        
        
        return N, tri

