import numpy as np
def normal(values, sigma, mu = 0):
    return -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (0 - values) ** 2 / sigma ** 2)