from dataclasses import dataclass
import numpy as np


@dataclass
class StructuredGrid:
    origin: np.ndarray
    maximum: np.ndarray
    step_vector: np.ndarray
    nsteps: np.ndarray
    data: np.ndarray

    def to_dict(self):
        return {
            "origin": self.origin,
            "maximum": self.maximum,
            "step_vector": self.step_vector,
            "nsteps": self.nsteps,
            "data": self.data,
        }

    def vtk(self):
        pass
