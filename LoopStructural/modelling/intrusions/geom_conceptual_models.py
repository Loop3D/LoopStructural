# Geometrical conceptual models for lateral and vertical extent of intrusions
import numpy as np
import pandas as pd
from ...utils import getLogger

logger = getLogger(__name__)


# conceptual geometrical models


def ellipse_function(
    lateral_contact_data=pd.DataFrame(),
    model=True,  # True to cover the extent of the model, regardless of data distribution
    minP=None,
    maxP=None,
    minS=None,
    maxS=None,
):

    if lateral_contact_data.empty:
        return model, minP, maxP, minS, maxS

    else:

        if minP is None:
            minP = lateral_contact_data["coord1"].min()
        if maxP is None:
            maxP = lateral_contact_data["coord1"].max()
        if minS is None:
            minS = lateral_contact_data["coord2"].abs().min()
        if maxS is None:
            maxS = lateral_contact_data["coord2"].max()

        a = (maxP - minP) / 2
        b = (maxS - minS) / 2

        po = minP + (maxP - minP) / 2

        p_locations = lateral_contact_data.loc[:, "coord1"].copy().to_numpy()

        s = np.zeros([len(p_locations), 2])

        s[np.logical_and(p_locations > minP, p_locations < maxP), 0] = b * np.sqrt(
            1
            - np.power(
                (p_locations[np.logical_and(p_locations > minP, p_locations < maxP)] - po) / a,
                2,
            )
        )
        s[np.logical_and(p_locations > minP, p_locations < maxP), 1] = -b * np.sqrt(
            1
            - np.power(
                (p_locations[np.logical_and(p_locations > minP, p_locations < maxP)] - po) / a,
                2,
            )
        )

        return s


def constant_function(
    othercontact_data=pd.DataFrame(),
    mean_growth=None,
    minP=None,
    maxP=None,
    minS=None,
    maxS=None,
    vertex=None,
):

    if othercontact_data.empty:
        return mean_growth

    if mean_growth is None:
        mean_growth = othercontact_data.loc[:, 'coord1'].mean()

    data_ps = np.array([othercontact_data.loc[:, 'coord1'], othercontact_data.loc[:, 'coord2']]).T

    conceptual_growth = np.ones([len(data_ps), 2]) * mean_growth

    return conceptual_growth


def obliquecone_function(
    othercontact_data=pd.DataFrame(),
    mean_growth=None,
    minP=None,
    maxP=None,
    minS=None,
    maxS=None,
    vertex=None,
):
    # import math

    if othercontact_data.empty:
        return mean_growth

    ps_locations = othercontact_data.loc[:, ['coord1', 'coord2']].to_numpy()

    minP = 1.5 * minP
    maxP = 1.5 * maxP
    minS = 1.5 * minS
    maxS = 1.5 * maxS

    a = (maxP - minP) / 2  # semi-major axis
    b = (maxS - minS) / 2  # semi-minor axis
    a2 = pow(a, 2)
    b2 = pow(b, 2)

    po = minP + a  # p coordinate of ellipsis centre
    so = minS + b  # s coordinate of ellipsis centre

    alpha = vertex[0]  # p coordinate of vertex
    beta = vertex[2]  # g coordinate of vertex
    gamma = vertex[1]  # l coordinate of vertex

    growth = np.zeros([len(ps_locations), 2])  # container for results

    p = ps_locations[:, 0]
    s = ps_locations[:, 1]

    A = alpha - po
    B = beta * (p[:] - alpha)
    C = gamma - so
    D = beta * (s[:] - gamma)

    F = pow(A * b, 2) + pow(C * a, 2) - a2 * b2
    G = 2 * (B * A * b2 + C * D * a2)
    H = pow(b * B, 2) + pow(a * D, 2)

    constant_g2 = F
    constant_g = -2 * F * beta - G
    constant_1 = F * pow(beta, 2) + G * beta + H

    discriminant = pow(constant_g, 2) - 4 * constant_g2 * constant_1
    discriminant[discriminant < 0] = 0

    growth[:, 0] = -(constant_g + np.sqrt(discriminant)) / (2 * constant_g2)
    growth[:, 1] = -(constant_g - np.sqrt(discriminant)) / (2 * constant_g2)

    return growth
