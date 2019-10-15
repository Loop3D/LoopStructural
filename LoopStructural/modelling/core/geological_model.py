from LoopStructural.modelling.core.geological_points import GPoint, IPoint
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator, GeologicalFeature
from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid

import numpy as np
import networkx as nx

import logging
logger = logging.getLogger(__name__)


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes the rescaling
    of the model between 0 and 1.
    """
    def __init__(self, origin, maximum):
        """

        Parameters
        ----------
        origin - numpy array specifying the origin of the model
        maximum - numpy array specifying the maximum extent of the model
        """
        self.graph = nx.DiGraph()
        self.features = {}
        self.data = {}


        # we want to rescale the model area so that the maximum length is
        # 1
        self.origin = origin
        self.maximum = maximum
        lengths = self.maximum - self.origin
        self.scale_factor = np.max(lengths)

        self.bounding_box = np.zeros((2, 3))
        self.bounding_box[0, :] = self.maximum-self.origin
        self.bounding_box /= self.scale_factor

    def get_interpolator(self, interpolatortype, nelements, buffer, region = "everywhere"):
        interpolator = None
        bb = np.copy(self.bounding_box)
        bb[0,:] -= buffer
        bb[1,:] += buffer
        if interpolatortype == "PLI":
            mesh = TetMesh()

            mesh.setup_mesh(bb, n_tetra=nelements,)
            return PLI(mesh)

        if interpolatortype == 'FDI':
            # number of elements should be divided roughly to match the shape
            ratio = bb[1,:] / np.sum(bb[1,:])
            ratio *= nelements
            ratio = ratio.astype(int)
            step_vector = 1. / np.max(ratio)
            grid = StructuredGrid(nsteps=ratio,step_vector=step_vector)
            return FDI(grid)

    def add_structural_frame(self, frame):
        # self.features[frame.name] = frame
        self.graph.add_node(frame, name=frame.name)

    def add_fold(self, fold):
        self.graph.add_node(fold, name=fold.name)

    def add_feature(self, feature, name):
        self.features[name] = feature

    def rescale(self, points):
        points*=self.scale_factor
        points+=self.origin
        return points

    def build_model(self):
        pass
