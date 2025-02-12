from abc import ABCMeta, abstractmethod
import numpy as np
from typing import Tuple


class BaseSupport(metaclass=ABCMeta):
    """
    Base support class
    """

    @abstractmethod
    def __init__(self):
        """
        This class is the base
        """

    @abstractmethod
    def evaluate_value(self, evaluation_points: np.ndarray, property_array: np.ndarray):
        """
        Evaluate the value of the support at the evaluation points
        """
        pass

    @abstractmethod
    def evaluate_gradient(self, evaluation_points: np.ndarray, property_array: np.ndarray):
        """
        Evaluate the gradient of the support at the evaluation points
        """
        pass

    @abstractmethod
    def inside(self, pos):
        """
        Check if a position is inside the support
        """
        pass

    @abstractmethod
    def onGeometryChange(self):
        """
        Called when the geometry changes
        """
        pass

    @abstractmethod
    def get_element_for_location(
        self, pos: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get the element for a location
        """
        pass

    @abstractmethod
    def get_element_gradient_for_location(
        self, pos: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        pass

    @property
    @abstractmethod
    def elements(self):
        """
        Return the elements
        """
        pass

    @property
    @abstractmethod
    def n_elements(self):
        """
        Return the number of elements
        """
        pass

    @property
    @abstractmethod
    def n_nodes(self):
        """
        Return the number of points
        """
        pass

    @property
    @abstractmethod
    def nodes(self):
        """
        Return the nodes
        """
        pass

    @property
    @abstractmethod
    def barycentre(self):
        """
        Return the number of dimensions
        """
        pass

    @property
    @abstractmethod
    def dimension(self):
        """
        Return the number of dimensions
        """
        pass

    @property
    @abstractmethod
    def element_size(self):
        """
        Return the element size
        """
        pass

    @abstractmethod
    def vtk(self, node_properties={}, cell_properties={}):
        """
        Return a vtk object
        """
        pass

    @abstractmethod
    def set_nelements(self, nelements) -> int:
        pass
