from FME.modelling.structural_frame import StructuralFrame
import numpy as np


class FoldFrame(StructuralFrame):
    """
    A structural frame that can calculate the fold axis/limb rotation angle
    """
    def __init__(self, name, features):
        super().__init__(name, features)

    def calculate_fold_axis_rotation(self, points):
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(s1g, axis=1)[:, None]
        # s1 = self.features[0].evaluate_value(points[:,:3])
        s1gyg = self.features[1].evaluate_gradient(points[:, :3])
        s1gyg /= np.linalg.norm(s1gyg, axis=1)[:, None]
        # s1gy = self.features[0].evaluate_value(points[:,:3])
        l1 = points[:,3:]#np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1, axis=1)[:, None]
        #einsum dot product
        far = np.einsum('ij,ij->i', s1gyg, l1)
        far = np.rad2deg(np.arccos(far))
        #scalar triple product
        stp = np.einsum('ij,ij->i', np.cross(l1, s1gyg, axisa=1, axisb=1), s1g)
        #check bounds
        # far[stp < 0] = 360.- far[stp < 0]
        # far[far>90] = far[far>90]+-180
        # far[far<-90] = far[far<-90]+180
        return far

    def calculate_fold_limb_rotation(self, points, axis):
        s0g = points[:,3:]
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(s1g,axis=1)[:,None]
        s1 = self.features[0].evaluate_value(points[:,:3])
        
        projected_s0 = s0g - np.einsum('ij,ij->i',axis,s0g)[:,None]*s0g
        projected_s1 = s1g - np.einsum('ij,ij->i',axis,s1g)[:,None]*s1g
        projected_s0/=np.linalg.norm(projected_s0, axis=1)[:,None]
        projected_s1/=np.linalg.norm(projected_s1, axis=1)[:,None]
        r2 = np.einsum('ij,ij->i', projected_s1, projected_s0)#s1g,s0g)#

        vv = np.cross(s1g, s0g, axisa=1, axisb=1)
        ds = np.einsum('ik,ij->i', axis, vv)
        flr = np.where(ds > 0, np.rad2deg(np.arcsin(r2)), (- np.rad2deg(np.arcsin(r2))))
        flr = np.where(flr < -90, (180.+flr), flr)
        flr = np.where(flr > 90, -(180.-flr), flr)
        return flr

    def calculate_intersection_lineation(self, points):
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(points[:,:3],axis=1)[:,None]
        s0g = points[:,3:]
        s0g /= np.linalg.norm(s0g,axis=1)[:,None]
        l1 = np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1,axis=1)[:,None]
        return l1
