# import scipy as sc
import scipy.stats as sct

import numpy as np
import pandas as pd

from ...utils import getLogger, rng

logger = getLogger(__name__)


def geometric_scaling_parameters(
    intrusion_type: str,
):
    """
    Get geometric scaling parameters for a given intrusion type

    Parameters
    ----------
    intrusion_type : str
        intrusion type

    Returns
    -------
    tuple(float, float, float, float)
        scaling parameters
    """
    geom_scaling_a_avg = {
        "plutons": 0.81,
        "laccoliths": 0.92,
        "major_mafic_sills": 0.85,
        "mesoscale_mafic_sills": 0.49,
        "minor_mafic_sills": 0.91,
    }
    geom_scaling_a_stdv = {
        "plutons": 0.12,
        "laccoliths": 0.11,
        "major_mafic_sills": 0.1,
        "mesoscale_mafic_sills": 0.13,
        "minor_mafic_sills": 0.25,
    }
    geom_scaling_b_avg = {
        "plutons": 1.08,
        "laccoliths": 0.12,
        "major_mafic_sills": 0.01,
        "mesoscale_mafic_sills": 0.47,
        "minor_mafic_sills": 0.27,
    }
    geom_scaling_b_stdv = {
        "plutons": 1.38,
        "laccoliths": 0.02,
        "major_mafic_sills": 0.02,
        "mesoscale_mafic_sills": 0.33,
        "minor_mafic_sills": 0.04,
    }

    a_avg = geom_scaling_a_avg.get(intrusion_type)
    a_stdv = geom_scaling_a_stdv.get(intrusion_type)
    b_avg = geom_scaling_b_avg.get(intrusion_type)
    b_stdv = geom_scaling_b_stdv.get(intrusion_type)

    return a_avg, a_stdv, b_avg, b_stdv


def thickness_from_geometric_scaling(length: float, intrusion_type: str) -> float:
    """Calculate thickness of intrusion using geometric scaling parameters

    Parameters
    ----------
    length : float
        intrusion length
    intrusion_type : str
        type of intrusion

    Returns
    -------
    float
        thickness of intrusion
    """

    a_avg, a_stdv, b_avg, b_stdv = geometric_scaling_parameters(intrusion_type)

    n_realizations = 10000
    maxT = 0
    a = sct.norm.ppf(rng.random(n_realizations), loc=a_avg, scale=a_stdv)
    b = sct.norm.ppf(rng.random(n_realizations), loc=b_avg, scale=b_stdv)
    maxT = b * np.power(length, a)
    maxT[maxT < 0] = None
    mean_t = np.nanmean(maxT)

    logger.info("Building intrusion of thickness {}".format(mean_t))

    return mean_t


def contact_pts_using_geometric_scaling(
    thickness: float, points_df: pd.DataFrame, inflation_vector: np.ndarray
):
    """Generate contact points for an intrusion using geometric scaling parameter and the
    inflation vector to translate the points

    Parameters
    ----------
    thickness : float
        intrusion thickness
    points_df : pd.DataFrame
        dataframe of contact points
    inflation_vector : np.ndarray
        inflation direction of the intrusion

    Returns
    -------
    tuple
        contact points
    """
    translation_vector = (
        inflation_vector
        / np.linalg.norm(inflation_vector, axis=1).reshape(1, len(inflation_vector)).T
    ) * thickness
    points_translated = points_df.loc[:, ["X", "Y", "Z"]].copy() + translation_vector
    points_translated_xyz = points_translated.to_numpy()

    return points_translated, points_translated_xyz
