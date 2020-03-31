import sys

from LoopStructural.utils.helper import get_vectors
from .geological_interpolator import GeologicalInterpolator

sys.path.append('/home/lgrose/dev/cpp/surfe/')
import surfepy
import logging
import numpy as np

logger = logging.getLogger(__name__)


class SurfeRBFInterpolator(GeologicalInterpolator):
    def __init__(self, method='horizons'):
        GeologicalInterpolator.__init__(self)
        self.surfe = None
        if method == 'single_surface':
            self.surfe = surfepy.Surfe_API(1)
        if method == 'Laujaunie':
            self.surfe = surfepy.Surfe_API(2)
        if method == 'horizons':
            logger.info("Using surfe horizon")
            self.surfe = surfepy.Surfe_API(4)

    def add_gradient_ctr_pts(self):
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            strike_vector, dip_vector = get_vectors(points[:, 3:6])

            strike_vector = np.hstack([points[:, :3], strike_vector.T])
            dip_vector = np.hstack([points[:, :3], dip_vector.T])
            self.surfe.SetTangentConstraints(strike_vector)
            self.surfe.SetTangentConstraints(dip_vector)
        pass

    def add_norm_ctr_pts(self):
        points = self.get_norm_constraints()
        if points.shape[0]>0:
            self.surfe.SetPlanarConstraints(points[:, :6])

    def add_ctr_pts(self):

        points = self.get_value_constraints()
        if points.shape[0]> 0:
            self.surfe.SetInterfaceConstraints(points[:,:4])
        # for i in range(points.shape[0]):
        #
        #     self.surfe.AddInterfaceConstraint(points[i, 0], points[i, 1], points[i, 2], points[i, 3], )

    def add_tangent_ctr_pts(self):
        pass

    def _solve(self, **kwargs):
        self.surfe.ComputeInterpolant()

    def _setup_interpolator(self, **kwargs):
        self.add_gradient_ctr_pts()
        self.add_norm_ctr_pts()
        self.add_ctr_pts()
        self.add_tangent_ctr_pts()

        kernel = kwargs.get("kernel", 'r3')
        logger.info("Setting surfe RBF kernel to %s" % kernel)
        regression = kwargs.get("regression_smoothing", 0.)
        if regression > 0:
            logger.info("Using regression smoothing %f" % regression)

        greedy = kwargs.get("greedy", (0, 0))

        if greedy[0] > 0 or greedy[1] > 0:
            logger.info("Using greedy algorithm: inferface %f and angular %f" %
                        (greedy[0], greedy[1]))
            self.surfe.SetRegressionSmoothin(True, greedy[0], greedy[1])

        global_anisotropy = kwargs.get("anisotropy", False)
        if global_anisotropy:
            logger.info("Using global anisotropy")
            self.surfe.SetGlobalAnisotropy(global_anisotropy)
    def update(self):
        return self.surfe.InterpolantComputed()