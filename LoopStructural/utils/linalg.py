import numpy as np


def normalise(v):
    v = np.array(v)

    np.linalg.norm(v, axis=1)
    return v / np.linalg.norm(v)
