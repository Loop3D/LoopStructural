"""
Piecewise linear interpolator
"""
import logging

import numpy as np

from LoopStructural.interpolators.discrete_interpolator import \
    DiscreteInterpolator
from LoopStructural.utils.helper import get_vectors

logger = logging.getLogger(__name__)


class P1Interpolator(DiscreteInterpolator):
    """

    """
    def __init__(self, mesh):
        """
        Piecewise Linear Interpolator
        Approximates scalar field by finding coefficients to a piecewise linear
        equation on a tetrahedral mesh. Uses constant gradient regularisation.

        Parameters
        ----------
        mesh - TetMesh
            interpolation support
        """

        self.shape = 'rectangular'
        DiscreteInterpolator.__init__(self, mesh)
        # whether to assemble a rectangular matrix or a square matrix
        self.interpolator_type = 'P1'
        self.nx = len(self.support.nodes[self.region])
        self.support = mesh

        self.interpolation_weights = {'cgw': 0.1, 'cpw': 1., 'npw': 1.,
                                      'gpw': 1., 'tpw': 1., 'ipw': 1.}
    
    def add_gradient_ctr_pts(self,w=1.0):
        pass
    
    def add_norm_ctr_pts(self,w=1.0):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:,:2])
            inside = elements  > -1
            area = self.support.element_area(elements[inside])
            wt = np.ones(area.shape[0])
            wt*=w*area
            # print(grad[inside,:,:].shape)
            # print(self.support.elements[elements[inside]].shape)
            elements = np.tile(self.support.elements[elements[inside]], (2, 1, 1))

            elements = elements.swapaxes(0,1)
            # elements = elements.swapaxes(0,2)
            grad = grad.swapaxes(1,2)
            
            self.add_constraints_to_least_squares(
                                                grad[inside,:,:]*wt[:,None,None],
                                                points[inside,3:5]*wt[:,None],
                                                elements,
                                                name = 'norm')

        pass

    def add_ctr_pts(self, w=1.0):
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            N, tri = self.support.evaluate_shape(points[:,:2])
            mask = tri > 0
            area = self.support.element_area(tri[mask])
            wt = np.ones(area.shape[0])
            wt*=w*area
            self.add_constraints_to_least_squares(N[mask,:]*wt[:,None],points[mask,3]*wt,
                                                  self.support.elements[tri[mask],:],
                                                  name='value')

    def minimize_edge_jumps(self,stren,w=0.1,maxmdDist=None): #NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        #iterate over all triangles
        # flag inidicate which triangles have had all their relationships added 
        v1 = self.support.nodes[self.support.edges][:,0,:]
        v2 = self.support.nodes[self.support.edges][:,1,:]
        ncp = 2
        # cp = np.zeros((v1.shape[0],ncp,2))
        # cp[:,0] = 0.25*v1 + 0.75*v2
        # cp[:,1] = 0.27*v1 + 0.25*v2
        bc_t1 = self.support.barycentre(self.support.edge_relationships[:,0])
        bc_t2 = self.support.barycentre(self.support.edge_relationships[:,1])

        v = v1-v2
        e_len = np.linalg.norm(v,axis=1)
        normal = np.array([v[:,1],-v[:,0]]).T
        normal /= np.linalg.norm(normal,axis=1)[:,None]
        # evaluate the shape function for the edges for each neighbouring triangle
        Dt, tri1 = self.support.evaluate_shape_derivatives(bc_t1,elements=self.support.edge_relationships[:,0])
        Dn, tri2 = self.support.evaluate_shape_derivatives(bc_t2,elements=self.support.edge_relationships[:,1])

        # constraint for each cp is triangle - neighbour create a Nx12 matrix 
        const_t = np.einsum('ij,ikj->ik',normal,Dt)
        const_n = -np.einsum('ij,ikj->ik',normal,Dn)
        # const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        # const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)

        const = np.hstack([const_t,const_n])
        
        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[tri1],self.support.elements[tri2]])
        # tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[tri2]])
        # add cp1 and cp2 to the least squares system
        self.add_constraints_to_least_squares(const*e_len[:,None]*w,np.zeros(const.shape[0]),tri_cp1, name='edge jump cp1')
        # p2.add_constraints_to_least_squares(const_cp2*e_len[:,None]*w,np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')
    
    def minimize_gradient_norm(self,stren,w=0.1,maxmdDist=None): #NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        #iterate over all triangles
        # flag inidicate which triangles have had all their relationships added 
        v1 = self.support.nodes[self.support.edges][:,0,:]
        v2 = self.support.nodes[self.support.edges][:,1,:]
        ncp = 2
        # cp = np.zeros((v1.shape[0],ncp,2))
        # cp[:,0] = 0.25*v1 + 0.75*v2
        # cp[:,1] = 0.27*v1 + 0.25*v2
        bc_t1 = self.support.barycentre(self.support.edge_relationships[:,0])
        bc_t2 = self.support.barycentre(self.support.edge_relationships[:,1])

        v = v1-v2
        e_len = np.linalg.norm(v,axis=1)
        normal = np.array([v[:,1],-v[:,0]]).T
        normal /= np.linalg.norm(normal,axis=1)[:,None]
        # evaluate the shape function for the edges for each neighbouring triangle
        Dt, tri1 = self.support.evaluate_shape_derivatives(bc_t1,elements=self.support.edge_relationships[:,0])
        Dn, tri2 = self.support.evaluate_shape_derivatives(bc_t2,elements=self.support.edge_relationships[:,1])

        # constraint for each cp is triangle - neighbour create a Nx12 matrix 
      




        # const_t = np.einsum('ij,ikj->ik',normal,Dt)
        # const_n = -np.einsum('ij,ikj->ik',normal,Dn)
        # const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        # const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)

        const = np.hstack([Dt,-Dn]).T
        const = const.swapaxes(1,2)

        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[tri1],self.support.elements[tri2]])
        tri_cp1 = np.tile(tri_cp1,(2,1,1))
        B = np.zeros((const.shape[0],const.shape[1]))
        # np.tile(tri_cp1,(2,1,1))
        print(tri_cp1.shape,const.shape,B.shape)
        # tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[tri2]])
        # add cp1 and cp2 to the least squares system
        self.add_constraints_to_least_squares(const*e_len[None,:,None]*w,B,tri_cp1, name='constant norm')
        # p2.add_constraints_to_least_squares(const_cp2*e_len[:,None]*w,np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')
        

            