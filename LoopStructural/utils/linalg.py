import numpy as np


def normalise(v):
    v = np.array(v)

    norm = np.linalg.norm(v, axis=1)
    return v / np.linalg.norm(v)
