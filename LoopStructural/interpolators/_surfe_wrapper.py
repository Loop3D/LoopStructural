"""
Wrapper for using surfepy
"""

from ..utils.helper import get_vectors
from ..interpolators import GeologicalInterpolator

import numpy as np

from ..utils import getLogger

logger = getLogger(__name__)
import surfepy


class SurfeRBFInterpolator(GeologicalInterpolator):
    """ """

    def __init__(self, method="single_surface"):
        GeologicalInterpolator.__init__(self)
        self.surfe = None
        if method == "single_surface":
            logger.info("Using single surface interpolator")
            self.surfe = surfepy.Surfe_API(1)
        if method == "Lajaunie" or method == "increments":
            logger.info("Using Lajaunie method")
            self.surfe = surfepy.Surfe_API(2)
        if method == "horizons":
            logger.info("Using surfe horizon")
            self.surfe = surfepy.Surfe_API(4)

    def add_gradient_ctr_pts(self):
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            logger.info("Adding ")
            strike_vector, dip_vector = get_vectors(points[:, 3:6])

            strike_vector = np.hstack([points[:, :3], strike_vector.T])
            dip_vector = np.hstack([points[:, :3], dip_vector.T])
            self.surfe.SetTangentConstraints(strike_vector)
            self.surfe.SetTangentConstraints(dip_vector)

    def add_norm_ctr_pts(self):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            self.surfe.SetPlanarConstraints(points[:, :6])

    def add_ctr_pts(self):

        points = self.get_value_constraints()
        if points.shape[0] > 0:
            # self.surfe.SetInterfaceConstraints(points[:,:4])
            for i in range(points.shape[0]):

                self.surfe.AddInterfaceConstraint(
                    points[i, 0],
                    points[i, 1],
                    points[i, 2],
                    points[i, 3],
                )

    def add_tangent_ctr_pts(self):
        points = self.get_tangent_constraints()
        if points.shape[0] > 0:
            self.surfe.SetTangentConstraints(points[:, :6])

    def _solve(self, **kwargs):
        self.surfe.ComputeInterpolant()

    def _setup_interpolator(self, **kwargs):
        """
        Setup the interpolator

        Parameters
        ----------
        kernel: str
            kernel for interpolation r3, r, Gaussian, Multiquadrics, Inverse Multiquadrics
            Thin Plate Spline, WendlandC2, MaternC4
        regression: float
            smoothing parameter default 0
        greedy: tuple
            greedy parameters first is interface threshold, second is angular threshold
             default (0,0)
        poly_order: int
            order of the polynomial used for interpolation, default 1
        radius: float
            radius of the kernel, default None but required for SPD kernels
        anisotropy: bool
            apply global anisotropy from eigenvectors of orientation constraints, default False


        """
        self.add_gradient_ctr_pts()
        self.add_norm_ctr_pts()
        self.add_ctr_pts()
        self.add_tangent_ctr_pts()

        kernel = kwargs.get("kernel", "r3")
        logger.info("Setting surfe RBF kernel to %s" % kernel)
        self.surfe.SetRBFKernel(kernel)
        regression = kwargs.get("regression_smoothing", 0.0)
        if regression > 0:
            logger.info("Using regression smoothing %f" % regression)
            self.surfe.SetRegressionSmoothing(True, regression)
        greedy = kwargs.get("greedy", (0, 0))

        if greedy[0] > 0 or greedy[1] > 0:
            logger.info(
                "Using greedy algorithm: inferface %f and angular %f" % (greedy[0], greedy[1])
            )
            self.surfe.SetGreedyAlgorithm(True, greedy[0], greedy[1])
        poly_order = kwargs.get("poly_order", None)
        if poly_order:
            logger.info("Setting poly order to %i" % poly_order)
            self.surfe.SetPolynomialOrder(poly_order)
        global_anisotropy = kwargs.get("anisotropy", False)
        if global_anisotropy:
            logger.info("Using global anisotropy")
            self.surfe.SetGlobalAnisotropy(global_anisotropy)
        radius = kwargs.get("radius", False)
        if radius:
            logger.info("Setting RBF radius to %f" % radius)
            self.surfe.SetRBFShapeParameter(radius)

    def update(self):
        return self.surfe.InterpolantComputed()

    def evaluate_value(self, evaluation_points):
        """Evaluate surfe interpolant at points

        Parameters
        ----------
        evaluation_points : array of locations N,3
            xyz of locations to evaluate

        Returns
        -------
        np.array (N)
            value of interpolant at points
        """
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.surfe.EvaluateInterpolantAtPoints(evaluation_points[~mask])
        return evaluated

    def evaluate_gradient(self, evaluation_points):
        """Evaluate surfe interpolant gradient at points

        Parameters
        ----------
        evaluation_points : array of locations N,3
            xyz of locations to evaluate

        Returns
        -------
        np.array (N,3)
            gradient of interpolant at points

        """
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape)
        mask = np.any(evaluation_points == np.nan, axis=1)
        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask, :] = self.surfe.EvaluateVectorInterpolantAtPoints(
                evaluation_points[~mask]
            )
        return

    @property
    def nx(self):
        return self.get_data_locations().shape[0]
