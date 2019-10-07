from LoopStructural.modelling.core.geological_points import GPoint, IPoint
import numpy as np
import networkx as nx

import logging
logger = logging.getLogger(__name__)


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes the
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
        self.bounding_box[0, :] = self.maxmimum-self.origin
        self.bounding_box /= self.scale_factor

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
