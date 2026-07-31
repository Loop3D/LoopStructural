from abc import ABC, abstractmethod

import numpy as np


class BaseRepresentation(ABC):
    @abstractmethod
    def to_dict(self):
        pass

    @classmethod
    @abstractmethod
    def from_dict(cls, data):
        pass

    def __repr__(self):
        return f"{self.__class__.__name__}({self.to_dict()})"

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if not isinstance(other, BaseRepresentation):
            return NotImplemented
        return self.to_dict() == other.to_dict()

    def __ne__(self, other):
        eq_result = self.__eq__(other)
        if eq_result is NotImplemented:
            return NotImplemented
        return not eq_result

    def __hash__(self):
        return hash(tuple(sorted(self.to_dict().items())))

    @abstractmethod
    def evaluate_value(self, position: np.ndarray):
        raise NotImplementedError("Value evaluation not implemented for this representation")

    @abstractmethod
    def evaluate_gradient(self, position: np.ndarray):
        raise NotImplementedError("Gradient evaluation not implemented for this representation")

    def surfaces(self, value):
        raise NotImplementedError("Surface extraction not implemented for this representation")
