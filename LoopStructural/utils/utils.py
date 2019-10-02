import numpy as np
##TODO add global valiable for debug
def fme_debug(value):
    print(value)


def strike_symbol(strike):
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))

    vec = np.array([0, 1])
    rotated = R @ vec
    vec2 = np.array([-0.5, 0])
    r2 = R @ vec2
    return rotated, r2