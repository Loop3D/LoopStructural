import logging

import numpy as np

logger = logging.getLogger(__name__)


class FoldEvent:
    """

    """
    def __init__(self, foldframe, fold_axis_rotation=None, fold_limb_rotation=None, fold_axis=None, name='Fold'):
        """
        A fold event describes the geometry of the fold using a fold frame
        and two geometrical descriptors

        Parameters
        ----------
        foldframe the fold frame object
        fold_axis_rotation function for the fold axis rotation angle
        fold_limb_rotation function for the fold limb rotation angle
        """

        self.foldframe = foldframe
        self.fold_axis_rotation = fold_axis_rotation
        self.fold_limb_rotation = fold_limb_rotation
        self.fold_axis = fold_axis
        self.name = name

    def get_fold_axis_orientation(self, points):
        """
        gets the fold axis orientation for evaluation points

        Parameters
        ----------
        points locations to calculate fold axis numpy array Nx3

        Returns
        -------

        """
        # evaluate fold axis from rotation angle
        if self.fold_axis_rotation is not None:
            logger.info("Using fold_axis_rotation function")
            # get the gz direction
            dgx = self.foldframe.features[0].evaluate_gradient(points)
            dgy = self.foldframe.features[1].evaluate_gradient(points)
            dgx /= np.linalg.norm(dgx, axis=1)[:, None]
            dgy /= np.linalg.norm(dgy, axis=1)[:, None]
            # get gy
            gy = self.foldframe.features[1].evaluate_value(points)
            R1 = self.rot_mat(-dgx, self.fold_axis_rotation(gy))
            fold_axis = np.einsum('ijk,ki->kj', R1, dgy)
            fold_axis /= np.linalg.norm(fold_axis, axis=1)[:, None]
            return fold_axis

        # use constant fold axis
        if self.fold_axis is not None:
            logger.info("Using constant fold axis")
            return np.tile(self.fold_axis, (points.shape[0], 1))

    def get_deformed_orientation(self, points):
        """
        Calculate the normal to the folded foliation at locations

        Parameters
        ----------
        points - np.array
            location Nx3 array of x,y,z locations to evaluate fold

        Returns
        -------

        """
        fold_axis = self.get_fold_axis_orientation(points)
        gx = self.foldframe.features[0].evaluate_value(points)
        dgx = self.foldframe.features[0].evaluate_gradient(points)
        dgx /= np.linalg.norm(dgx, axis=1)[:, None]
        dgz = np.cross(dgx,fold_axis,axisa=1,axisb=1)
        # dgz = self.foldframe.features[2].evaluate_gradient(points)
        dgz /= np.linalg.norm(dgz, axis=1)[:, None]

        R2 = self.rot_mat(fold_axis, self.fold_limb_rotation(gx))
        fold_direction = np.einsum('ijk,ki->kj', R2, dgx)
        fold_direction /= np.sum(fold_direction, axis=1)[:, None]
        # calculate dot product between fold_direction and axis
        # if its less than 0 then inverse dgz
        d = np.einsum('ij,ik->i', fold_direction, fold_axis)
        dgz[d < 0] = -dgz[d < 0]
        return fold_direction, fold_axis, dgz

    def get_regularisation_direction(self, points):
        self.foldframe.features[2].evaluate_gradient(points)

    def rot_mat(self, axis, angle):
        """
        Create a rotation matrix for axis and angle

        Parameters
        ----------
        axis Nx3 vector for axis
        angle N array for angle in degrees

        Returns 3,3,N rotation matrix
        -------

        """
        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))
        C = 1.0 - c
        x = axis[:, 0]
        y = axis[:, 1]
        z = axis[:, 2]
        xs = x * s
        ys = y * s
        zs = z * s
        xC = x * C
        yC = y * C
        zC = z * C
        xyC = x * yC
        yzC = y * zC
        zxC = z * xC
        rotation_mat = np.zeros((3, 3, len(angle)))
        rotation_mat[0, 0, :] = x * xC + c
        rotation_mat[0, 1, :] = xyC - zs
        rotation_mat[0, 2, :] = zxC + ys

        rotation_mat[1, 0, :] = xyC + zs
        rotation_mat[1, 1, :] = y * yC + c
        rotation_mat[1, 2, :] = yzC - xs

        rotation_mat[2, 0, :] = zxC - ys
        rotation_mat[2, 1, :] = yzC + xs
        rotation_mat[2, 2, :] = z * zC + c
        return rotation_mat
