"""
Piecewise linear interpolator
"""
import logging

import numpy as np

from LoopStructural.interpolators.discrete_interpolator import \
    DiscreteInterpolator
from LoopStructural.utils.helper import get_vectors

logger = logging.getLogger(__name__)


class P2Interpolator(DiscreteInterpolator):
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
        self.interpolator_type = 'P2'
        self.nx = len(self.support.nodes[self.region])
        self.support = mesh

        self.interpolation_weights = {'cgw': 0.1, 'cpw': 1., 'npw': 1.,
                                      'gpw': 1., 'tpw': 1., 'ipw': 1.}
    
    def add_gradient_ctr_pts(self,w=1.0):
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:,:2])
            inside = elements  > -1
            area = self.support.element_area(elements[inside])
            wt = np.ones(area.shape[0])
            wt*=w*area
            A = np.einsum('ikj,ij->ik', grad[:,:], points[:,3:5])
            #A = np.einsum('ij,ikj->ik', grad[:,:], points[:,3:5])
            A *= area[idx, None]
# # A *= fold_orientation
            B = np.zeros(A.shape[0])
            idc = p2_mesh.elements[idx,:]
            self.add_constraints_to_least_squares(A, B, elements)
# p2.add_constraints_to_least_squares(A, B, idc,'fold ori')
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
    
    def add_gradient_orthogonal_constraint(self, points, vector, w=1.0,
                                           B=0):
        """
        constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        position
        normals
        w
        B

        Returns
        -------

        """
        pass
        # if points.shape[0] > 0:
        #     vertices, element_gradients, tetras, inside = self.support.get_tetra_gradient_for_location(points[:,:3])
        #     #e, inside = self.support.elements_for_array(points[:, :3])
        #     #nodes = self.support.nodes[self.support.elements[e]]
        #     vector /= np.linalg.norm(vector,axis=1)[:,None]
        #     vecs = vertices[:, 1:, :] - vertices[:, 0, None, :]
        #     vol = np.abs(np.linalg.det(vecs))  # / 6
        #     # d_t = self.support.get_elements_gradients(e)
        #     norm = np.linalg.norm(element_gradients, axis=2)
        #     element_gradients /= norm[:, :, None]

        #     A = np.einsum('ij,ijk->ik', vector, element_gradients)

        #     A *= vol[:, None]

        #     gi = np.zeros(self.support.n_nodes).astype(int)
        #     gi[:] = -1
        #     gi[self.region] = np.arange(0, self.nx).astype(int)
        #     w /= 3
        #     idc = gi[tetras]
        #     B = np.zeros(idc.shape[0])+B
        #     outside = ~np.any(idc == -1, axis=1)
        #     self.add_constraints_to_least_squares(A[outside, :] * w,
        #                                           B[outside], idc[outside, :])
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
    

    
    def minimise_grad_steepness(self, stren =0.0, w=0.1,maskall=False,wtfunc=None):
        """This constraint minimises the second derivative of the gradient
        mimimising the 2nd derivative should prevent high curvature solutions
        It is not added on the borders 

        Parameters
        ----------
        sten : float, optional
            [description], by default 0.0
        w : float, optional
            [description], by default 0.1
        """
        tri = np.arange(0,len(self.support.elements))
        mask = np.ones(self.support.neighbours.shape[0],dtype=bool)
        if maskall == False:
            mask[:] = np.all(self.support.neighbours > 0,axis=1)
        tri_points = self.support.nodes[self.support.elements[tri[mask],:],:]
        barycentre = np.mean(tri_points,axis=1)
        M_t = np.ones((tri_points.shape[0],3,3))
        M_t[:,:,1:] = tri_points[:,:3,:]
        area = np.abs(np.linalg.det(M_t))*0.5
        xyConst=self.support.evaluate_mixed_derivative(tri[mask])
        xxConst,yyConst=self.support.evaluate_shape_d2(tri[mask])
        wt = np.ones(tri_points.shape[0])
        wt  *= w*area

        d = np.linalg.norm((self.get_data_locations()[:,None,:2]-barycentre[None,:,:]),axis=2)
        min_dist = np.min(d,axis=0)
        min_dist/=np.max(min_dist)      
        wt*=(1+stren*min_dist)*area
        if wtfunc:
            wt=wtfunc(barycentre)*area
        idc = self.support.elements[tri[mask]]
        self.add_constraints_to_least_squares(xyConst*4*wt[:,None],np.zeros(xyConst.shape[0]),idc)
        
        self.add_constraints_to_least_squares(xxConst*wt[:,None],np.zeros(xxConst.shape[0]),idc)
        self.add_constraints_to_least_squares(yyConst*wt[:,None],np.zeros(yyConst.shape[0]),idc)

            
        
    def minimize_edge_jumps(self,stren,w=0.1,maxmdDist=None,wtfunc=None, vector_func=None): #NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        #iterate over all triangles
        # flag inidicate which triangles have had all their relationships added 
        v1 = self.support.nodes[self.support.edges][:,0,:]
        v2 = self.support.nodes[self.support.edges][:,1,:]
        ncp = 2
        cp = np.zeros((v1.shape[0],ncp,2))
        cp[:,0] = 0.25*v1 + 0.75*v2
        cp[:,1] = 0.27*v1 + 0.25*v2

        d = np.linalg.norm((self.get_data_locations()[:,None,:2]-cp[None,:,0,:]),axis=2)
        cp1_min_dist = np.min(d,axis=0)

        d = np.linalg.norm((self.get_data_locations()[:,None,:2]-cp[None,:,1,:]),axis=2)
        cp2_min_dist = np.min(d,axis=0)
        
        cp1_min_dist/=np.max([cp1_min_dist,cp2_min_dist]) 
        cp2_min_dist/=np.max([cp1_min_dist,cp2_min_dist]) 

        v = v1-v2
        e_len = np.linalg.norm(v,axis=1)
        normal = np.array([v[:,1],-v[:,0]]).T
        normal /= np.linalg.norm(normal,axis=1)[:,None]
        
        # evaluate normal if using vector func for cp1 
        if vector_func:
            normal = vector_func(cp[:,0])
        
        # evaluate the shape function for the edges for each neighbouring triangle
        cp1_Dt, cp1_tri1 = self.support.evaluate_shape_derivatives(cp[:,0],elements=self.support.edge_relationships[:,0])
        cp1_Dn, cp1_tri2 = self.support.evaluate_shape_derivatives(cp[:,0],elements=self.support.edge_relationships[:,1])
        
        # constraint for each cp is triangle - neighbour create a Nx12 matrix 
        const_t_cp1 = np.einsum('ij,ikj->ik',normal,cp1_Dt)
        const_n_cp1 = -np.einsum('ij,ikj->ik',normal,cp1_Dn)
        
        # evaluate normal if using vector func for cp2 
        if vector_func:
            normal = vector_func(cp[:,1])
        # evaluate the shape function for the edges for each neighbouring triangle
        cp2_Dt, cp2_tri1 = self.support.evaluate_shape_derivatives(cp[:,1],elements=self.support.edge_relationships[:,0])
        cp2_Dn, cp2_tri2 = self.support.evaluate_shape_derivatives(cp[:,1],elements=self.support.edge_relationships[:,1])
        # constraint for each cp is triangle - neighbour create a Nx12 matrix 

        const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)

        const_cp1 = np.hstack([const_t_cp1,const_n_cp1])
        const_cp2 = np.hstack([const_t_cp2,const_n_cp2])
        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[cp1_tri1],self.support.elements[cp1_tri2]])
        tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[cp2_tri2]])
        # add cp1 and cp2 to the least squares system
        wt = np.zeros(tri_cp1.shape[0])
        wt[:] = w
        if wtfunc:
            wt=wtfunc(tri_cp1)
        self.add_constraints_to_least_squares(const_cp1*e_len[:,None]*wt[:,None],np.zeros(const_cp1.shape[0]),tri_cp1, name='edge jump cp1')
        if wtfunc:
            wt=wtfunc(tri_cp2)
        self.add_constraints_to_least_squares(const_cp2*e_len[:,None]*wt[:,None],np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')

    def minimize_edge_jump_magnitude(self,stren,w=0.1,maxmdDist=None,wtfunc=None, vector_func=None): #NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        #iterate over all triangles
        # flag inidicate which triangles have had all their relationships added 
        v1 = self.support.nodes[self.support.edges][:,0,:]
        v2 = self.support.nodes[self.support.edges][:,1,:]
        ncp = 2
        midedge = v1-v2
        
        cp = np.zeros((v1.shape[0],ncp,2))
        cp[:,0] = 0.25*v1 + 0.75*v2
        cp[:,1] = 0.27*v1 + 0.25*v2

        d = np.linalg.norm((self.get_data_locations()[:,None,:2]-cp[None,:,0,:]),axis=2)
        cp1_min_dist = np.min(d,axis=0)

        d = np.linalg.norm((self.get_data_locations()[:,None,:2]-cp[None,:,1,:]),axis=2)
        cp2_min_dist = np.min(d,axis=0)
        
        cp1_min_dist/=np.max([cp1_min_dist,cp2_min_dist]) 
        cp2_min_dist/=np.max([cp1_min_dist,cp2_min_dist]) 

        v = v1-v2
        e_len = np.linalg.norm(v,axis=1)
        normal = np.array([v[:,1],-v[:,0]]).T
        normal /= np.linalg.norm(normal,axis=1)[:,None]
        
        # evaluate normal if using vector func for cp1 
        if vector_func:
            normal = vector_func(cp[:,0])
        
        # evaluate the shape function for the edges for each neighbouring triangle
        cp1_Dt, cp1_tri1 = self.support.evaluate_shape_derivatives(cp[:,0],elements=self.support.edge_relationships[:,0])
        cp1_Dn, cp1_tri2 = self.support.evaluate_shape_derivatives(cp[:,0],elements=self.support.edge_relationships[:,1])
        
        # constraint for each cp is triangle - neighbour create a Nx12 matrix 
        const_t_cp1 = np.einsum('ij,ikj->ik',normal,cp1_Dt)
        const_n_cp1 = -np.einsum('ij,ikj->ik',normal,cp1_Dn)
        
        # evaluate normal if using vector func for cp2 
        if vector_func:
            normal = vector_func(cp[:,1])
        # evaluate the shape function for the edges for each neighbouring triangle
        cp2_Dt, cp2_tri1 = self.support.evaluate_shape_derivatives(cp[:,1],elements=self.support.edge_relationships[:,0])
        cp2_Dn, cp2_tri2 = self.support.evaluate_shape_derivatives(cp[:,1],elements=self.support.edge_relationships[:,1])
        # constraint for each cp is triangle - neighbour create a Nx12 matrix 

        const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)

        const_cp1 = np.hstack([const_t_cp1,const_n_cp1])
        const_cp2 = np.hstack([const_t_cp2,const_n_cp2])
        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[cp1_tri1],self.support.elements[cp1_tri2]])
        tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[cp2_tri2]])
        # add cp1 and cp2 to the least squares system
        wt = np.zeros(tri_cp1.shape[0])
        wt[:] = w
        if wtfunc:
            wt=wtfunc(tri_cp1)
        self.add_constraints_to_least_squares(const_cp1*e_len[:,None]*wt[:,None],np.zeros(const_cp1.shape[0]),tri_cp1, name='edge jump cp1')
        if wtfunc:
            wt=wtfunc(tri_cp2)
        self.add_constraints_to_least_squares(const_cp2*e_len[:,None]*wt[:,None],np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')

    def evaluate_d2(self, evaluation_points):
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_d2(
                evaluation_points[~mask], self.propertyname)
        return evaluated
            