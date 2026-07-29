"""
Finite difference masks
"""

import numpy as np

from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


class Operator(object):
    """
    Finite difference masks for adding constraints for the derivatives and second derivatives
    Operator.Dx_mask gives derivative in x direction
    """

    z = np.zeros((3, 3))
    Dx_mask = np.array([z, [[0.0, 0.0, 0.0], [-0.5, 0.0, 0.5], [0.0, 0.0, 0.0]], z])
    Dy_mask = Dx_mask.swapaxes(1, 2)
    Dz_mask = Dx_mask.swapaxes(0, 2)

    Dx_forward_mask = np.array([z, [[0.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 0.0]], z])
    Dx_backward_mask = np.array([z, [[0.0, 0.0, 0.0], [0.0, -1.0, 1.0], [0.0, 0.0, 0.0]], z])
    Dy_forward_mask = Dx_forward_mask.swapaxes(1, 2)
    Dy_backward_mask = Dx_backward_mask.swapaxes(1, 2)
    Dz_forward_mask = Dx_forward_mask.swapaxes(0, 2)
    Dz_backward_mask = Dx_backward_mask.swapaxes(0, 2)

    Dxx_mask = np.array([z, [[0, 0, 0], [1, -2, 1], [0, 0, 0]], z])
    Dyy_mask = Dxx_mask.swapaxes(1, 2)
    Dzz_mask = Dxx_mask.swapaxes(0, 2)

    Dxy_mask = np.array([z, [[-0.25, 0, 0.25], [0, 0, 0], [0.25, 0, -0.25]], z]) / np.sqrt(2)
    Dxz_mask = Dxy_mask.swapaxes(0, 1)
    Dyz_mask = Dxy_mask.swapaxes(0, 2)

    # from https://en.wikipedia.org/wiki/Discrete_Laplace_operator
    Lapacian = np.array(
        [
            [[0, 0, 0], [0, 1, 0], [0, 0, 0]],  # first plane
            [[0, 1, 0], [1, -6, 1], [0, 1, 0]],  # second plane
            [[0, 0, 0], [0, 1, 0], [0, 0, 0]],  # third plane
        ]
    )
