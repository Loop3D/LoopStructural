from .dsi_helper import fold_cg, cg_cstr_to_coo_sq, cg_cstr_to_coo_sq_B
import numpy as np
class FoldEvent():
    def __init__(self,foldframe,fold_axis_rotation,fold_limb_rotation):
        self.foldframe = foldframe
        self.fold_axis_rotation = fold_axis_rotation
        self.fold_limb_rotation = fold_limb_rotation
        
    def get_fold_axis_orientation(self,points):
        #get the gz direction
        dgz = self.foldframe.get_gz(points,grad=True)
        dgy = self.foldframe.get_gy(points,grad=True)
        #get gy
        gy = self.foldframe.get_gy(points,grad=False)
        R1 = self.rot_mat(dgz,self.fold_axis_rotation(gy))
        fold_axis = np.einsum('ijk,ki->kj',R1,dgy)
        fold_axis/=np.sum(fold_axis,axis=1)[:,None]
        return fold_axis
    def get_deformed_orientation(self,points):
        fold_axis=self.get_fold_axis_orientation(points)
        gx = self.foldframe.get_gx(points,grad=False)
        dgx = self.foldframe.get_gx(points,grad=True)
        dgz = self.foldframe.get_gz(points,grad=True)
        R2 = self.rot_mat(fold_axis,self.fold_limb_rotation(gx))
        R2R = np.einsum('ijk,ki->kj',R2,dgx)
        R2R/=np.sum(R2R,axis=1)[:,None]
        return R2R,fold_axis, dgz    
    def get_regulariation_direction(self,points):
        self.foldframe.get_gz(points,grad=True)
    def rot_mat(self,axis,angle):
        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))
        C = 1.0 - c
        x = axis[:,0]
        y = axis[:,1]
        z = axis[:,2]
        xs = x*s
        ys = y*s
        zs = z*s
        xC = x*C
        yC = y*C
        zC = z*C
        xyC = x*yC
        yzC = y*zC
        zxC = z*xC
        rotation_mat = np.zeros((3,3,len(angle)))
        rotation_mat[0,0,:] = x*xC+c
        rotation_mat[0,1,:] = xyC-zs
        rotation_mat[0,2,:] = zxC+ys

        rotation_mat[1,0,:] = xyC+zs
        rotation_mat[1,1,:] = y*yC+c
        rotation_mat[1,2,:] = yzC-xs

        rotation_mat[2,0,:] = zxC -ys
        rotation_mat[2,1,:] = yzC+xs
        rotation_mat[2,2,:] = z*zC+c
        return rotation_mat    
class DSIFoldConstraint():
    """
    Calculates the DSI constraints for the fold using a fold event and mesh.
    """
    def __init__(self,fold,mesh,shape='rectangular'):
        self.mesh = mesh
        self.fold = fold
        self.shape = shape
        self.deformed_orientation = False
        self.fold_regularization = False
        self.fold_normalisation = False
        self.fold_axis = False
    def active_constraints(self):
        """
        Print out which fold constraints are turned on and their weights
        """
        pass 
    def use_deformed_orientation_cnstrt(self,w=1.):
        self.deformed_orientation = True
        self.deformed_orientation_w = w
        return
    def use_fold_axis_cnstrt(self,w=1.):
        self.fold_axis = True
        self.fold_axis_w = w
        return
    def use_regularisation_cnstrt(self,w=1.):
        self.fold_regularization = True
        self.fold_regularization_w = w
        return
    def use_normalisation_cnstrt(self,norm=1.,w=1.):
        self.fold_normalisation = True
        self.fold_normalisation_w = w
        self.fold_apparent_norm = norm
        return
    def get_constraints(self):
        A = []
        B = []
        if self.shape == 'square':
            B = np.zeros(self.mesh.n_nodes)
        row = []
        col = []
        c_ = 0
        cstrs = {}
        eg = self.mesh.get_elements_gradients(np.arange(self.mesh.n_elements))
        deformed_orientation, fold_axis, dgz = self.fold.get_deformed_orientation(self.mesh.barycentre)
        if self.deformed_orientation:
            if self.shape == 'square':
                c = np.einsum('ij,ijk->ik',deformed_orientation,eg)
                a, r, c = cg_cstr_to_coo_sq(c,self.mesh.elements,self.mesh.n_elements)
                a = np.array(a)
                a*=self.deformed_orientation_w
                A.extend(a.tolist())
                row.extend(np.array(r).tolist())
                col.extend(np.array(c).tolist())
        if self.fold_axis:
            if self.shape == 'square':
                a, r, c = cg_cstr_to_coo_sq(np.einsum('ij,ijk->ik',fold_axis,eg),\
                                            self.mesh.elements,self.mesh.n_elements)
                a = np.array(a)
                a*=self.fold_axis_w
                A.extend(a.tolist())
                row.extend(np.array(r).tolist())
                col.extend(np.array(c).tolist())
        if self.fold_normalisation:
            if self.shape == 'square':
                Bb = np.zeros(self.mesh.n_elements)
                Bb[:] = self.fold_apparent_norm
                a, b, r, c = cg_cstr_to_coo_sq_B(np.einsum('ij,ijk->ik',dgz,eg), Bb,\
                                            self.mesh.elements,self.mesh.n_elements, self.mesh.n_nodes)
                
                a = np.array(a)
                a*=self.fold_normalisation_w
                A.extend(a.tolist())
                row.extend(np.array(r).tolist())
                col.extend(np.array(c).tolist())
                B+=b
        if self.fold_regularization:
            if self.shape == 'square':
                idc,c,ncons = fold_cg(eg,dgz,self.mesh.neighbours,self.mesh.elements,self.mesh.nodes)
                a,r,c = cg_cstr_to_coo_sq(c,idc,ncons)
                a = np.array(a)
                a*=self.fold_regularization_w
                A.extend(a.tolist())
                row.extend(np.array(r).tolist())
                col.extend(np.array(c).tolist())       
        return A,B,row,col
