from LoopStructural.modelling.features.structural_frame import StructuralFrame
import numpy as np

import logging
logger = logging.getLogger(__name__)


class FoldFrame(StructuralFrame):
    """
    A structural frame that can calculate the fold axis/limb rotation angle
    """
    def __init__(self, name, features):
        """
        Same constructor arguments as parent StructuralFrame
        Parameters
        ----------
        name
        features
        """
        super().__init__(name, features)

    def calculate_fold_axis_rotation(self, points):
        """
        Calculate the fold axis rotation angle by finding the angle between the
        intersection lineation and the gradient to the 1st coordinate of the fold frame
        Parameters
        ----------
        points

        Returns
        -------

        """
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(s1g, axis=1)[:, None]
        # s1 = self.features[0].evaluate_value(points[:,:3])
        s1gyg = self.features[1].evaluate_gradient(points[:, :3])
        s1gyg /= np.linalg.norm(s1gyg, axis=1)[:, None]
        # s1gy = self.features[0].evaluate_value(points[:,:3])
        l1 = points[:,3:]#np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1, axis=1)[:, None]

        # project s0 onto axis plane B X A X B
        projected_l1 = np.cross(s1g, np.cross(l1, s1g, axisa=1, axisb=1), axisa=1, axisb=1)
        projected_s1gyg = np.cross(s1g, np.cross(s1gyg, s1g, axisa=1, axisb=1), axisa=1, axisb=1)

        #einsum dot product
        far = np.einsum('ij,ij->i', projected_l1, projected_s1gyg)
        far = np.rad2deg(np.arccos(far))
        #scalar triple product
        stp = np.einsum('ij,ij->i', np.cross(l1, s1gyg, axisa=1, axisb=1), s1g)
        #check bounds
        # far[stp < 0] = 360.- far[stp < 0]
        # far[far>90] = far[far>90]+-180
        # far[far<-90] = far[far<-90]+180
        return far

    def calculate_fold_limb_rotation(self, points, axis):
        """
        Calculate the fold limb rotation angle using the axis specified and the normals to the folded foliation
        Parameters
        ----------
        points
        axis

        Returns
        -------

        """

        # get the normals from the points array
        s0g = points[:,3:]

        # calculate the gradient and value of the first coordinate of the fold frame
        # for the locations and normalise
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(s1g,axis=1)[:,None]
        s1 = self.features[0].evaluate_value(points[:,:3])

        # project s0 onto axis plane B X A X B
        projected_s0 = np.cross(axis, np.cross(s0g, axis, axisa=1, axisb=1), axisa=1, axisb=1)
        projected_s1 = np.cross(axis, np.cross(s1g, axis, axisa=1, axisb=1), axisa=1, axisb=1)
        projected_s0/=np.linalg.norm(projected_s0, axis=1)[:,None]
        projected_s1/=np.linalg.norm(projected_s1, axis=1)[:,None]
        r2 = np.einsum('ij,ij->i', s1g,s0g)#projected_s1, projected_s0)#
        # adjust the fold rotation angle so that its always between -90 and 90
        vv = np.cross(s1g, s0g, axisa=1, axisb=1)
        ds = np.einsum('ij,ij->i', axis, vv)
        flr = np.where(ds > 0, np.rad2deg(np.arcsin(r2)), (- np.rad2deg(np.arcsin(r2))))
        flr = np.where(flr < -90, (180.+flr), flr)
        flr = np.where(flr > 90, -(180.-flr), flr)
        return np.rad2deg(np.arcsin(r2))#flr

    def calculate_intersection_lineation(self, points):
        """
        Calculate the intersection lineation by finding the cross product between the first fold frame
        coordinate and the vector representing the normal to the folded foliation
        Parameters
        ----------
        points Nx6 array with x,y,z vx,vy,vz

        Returns Nx3 array of doubles
        -------

        """
        s1g = self.features[0].evaluate_gradient(points[:,:3])
        s1g /= np.linalg.norm(points[:,:3],axis=1)[:,None]
        s0g = points[:,3:]
        s0g /= np.linalg.norm(s0g,axis=1)[:,None]
        l1 = np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1,axis=1)[:,None]
        return l1
