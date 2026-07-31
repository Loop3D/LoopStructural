from abc import ABCMeta, abstractmethod
from typing import Tuple

import numpy as np


class BaseSupport(metaclass=ABCMeta):
    """
    Base support class
    """

    @abstractmethod
    def __init__(self):
        """
        This class is the base
        """

    def is_valid(self) -> bool:
        """
        Check if the support is valid
        """
        return True

    @abstractmethod
    def evaluate_value(self, evaluation_points: np.ndarray, property_array: np.ndarray):
        """
        Evaluate the value of the support at the evaluation points
        """

    @abstractmethod
    def evaluate_gradient(self, evaluation_points: np.ndarray, property_array: np.ndarray):
        """
        Evaluate the gradient of the support at the evaluation points
        """

    @abstractmethod
    def inside(self, pos):
        """
        Check if a position is inside the support
        """

    @abstractmethod
    def onGeometryChange(self):
        """
        Called when the geometry changes
        """

    @abstractmethod
    def get_element_for_location(
        self, pos: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get the element for a location
        """

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

    @property
    @abstractmethod
    def n_elements(self):
        """
        Return the number of elements
        """

    @property
    @abstractmethod
    def n_nodes(self):
        """
        Return the number of points
        """

    @property
    @abstractmethod
    def nodes(self):
        """
        Return the nodes
        """

    @property
    @abstractmethod
    def barycentre(self):
        """
        Return the number of dimensions
        """

    @property
    @abstractmethod
    def dimension(self):
        """
        Return the number of dimensions
        """

    @property
    @abstractmethod
    def element_size(self):
        """
        Return the element size
        """

    @abstractmethod
    def vtk(self, node_properties=None, cell_properties=None):
        """
        Return a vtk object
        """

    @abstractmethod
    def set_nelements(self, nelements) -> int:
        pass
