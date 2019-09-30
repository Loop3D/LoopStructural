import numpy as np

class Operator(object):

    z = np.zeros((3,3))
    Dx_mask = np.array([z,[
        [0.0,  0.0, 0.0],
        [-0.5, 0.0, 0.5],
        [0.0,  0.0, 0.0]],z
    ])
    Dy_mask = Dx_mask.swapaxes(1,2)
    Dz_mask = Dx_mask.swapaxes(0,2)

    Dxx_mask  = np.array([z,[
        [0, 0, 0],
        [1, -2, 1],
        [0, 0, 0]],z])
    Dyy_mask = Dxx_mask.swapaxes(1,2)
    Dzz_mask = Dxx_mask.swapaxes(0,2)

    Dxy_mask = np.array([z,[
        [-0.25, 0, 0.25],
        [0, 0, 0],
        [0.25, 0, -0.25]
    ],z])
    Dxz_mask = Dxy_mask.swapaxes(0,1)
    Dyz_mask = Dxy_mask.swapaxes(0,2)