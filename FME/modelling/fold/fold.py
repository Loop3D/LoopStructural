import numpy as np


class FoldEvent:
    """
    A fold event describes the geometry of the fold using a fold frame
    and two geometrical descriptors
    """
    def __init__(self, foldframe, fold_axis_rotation, fold_limb_rotation):
        """

        :param foldframe: the fold frame object
        :param fold_axis_rotation: function for the fold axis rotation angle
        :param fold_limb_rotation: function for the fold limb rotation angle
        """
        self.foldframe = foldframe
        self.fold_axis_rotation = fold_axis_rotation
        self.fold_limb_rotation = fold_limb_rotation
        
    def get_fold_axis_orientation(self, points):
        """
        gets the fold axis orientation for evaluation points
        :param points:
        :return:
        """

        # get the gz direction
        dgx = self.foldframe.features[0].evaluate_gradient(points)
        dgy = self.foldframe.features[1].evaluate_gradient(points)
        # get gy
        gy = self.foldframe.features[0].evaluate_value(points)
        R1 = self.rot_mat(dgx,self.fold_axis_rotation(gy))
        fold_axis = np.einsum('ijk,ki->kj',R1,dgy)
        fold_axis/=np.sum(fold_axis,axis=1)[:,None]
        return fold_axis

    def get_deformed_orientation(self, points):
        """
        Get the
        :param points:
        :return:
        """
        fold_axis=self.get_fold_axis_orientation(points)
        gx = self.foldframe.features[0].evaluate_value(points)
        dgx = self.foldframe.features[0].evaluate_gradient(points)
        dgz = self.foldframe.features[2].evaluate_gradient(points)
        R2 = self.rot_mat(fold_axis,self.fold_limb_rotation(gx))
        R2R = np.einsum('ijk,ki->kj',R2,dgx)
        R2R/=np.sum(R2R,axis=1)[:,None]
        return R2R,fold_axis, dgz

    def get_regulariation_direction(self, points):
        self.foldframe.features[2].evaluate_gradient(points)

    def rot_mat(self, axis, angle):
        """
        creates a rotation matrix from axis and angle
        :param axis:
        :param angle:
        :return:
        """
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
