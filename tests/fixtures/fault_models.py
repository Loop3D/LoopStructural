"""Small faulted models in the unit cube, for the fault block and kinematics tests.

The stratigraphic data are made the way LoopStructural evaluates a faulted
feature: each random present-day point is restored through the faults, and
the analytic field (z) of the restored point is its value.
"""

import numpy as np
import pandas as pd

from LoopStructural.modelling.features.fault._fault_function import (
    Composite,
    FaultDisplacement,
    Ones,
    Zeros,
)

COLUMNS = ["X", "Y", "Z", "val", "nx", "ny", "nz", "coord", "feature_name"]


def dip_vectors(dip_deg, dip_azimuth_deg):
    """Unit normal (up), down-dip and strike vectors of a plane."""
    d, a = np.radians(dip_deg), np.radians(dip_azimuth_deg)
    horizontal = np.array([np.sin(a), np.cos(a), 0.0])
    down_dip = np.cos(d) * horizontal + np.array([0.0, 0.0, -np.sin(d)])
    normal = np.sin(d) * horizontal + np.array([0.0, 0.0, np.cos(d)])
    return normal, down_dip, np.cross(normal, down_dip)


def constant_profile():
    """The hanging wall moves rigidly, the footwall does not move."""
    return FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)


def plane_points(center, dip, dip_azimuth, n=11, half=0.9):
    normal, down_dip, strike = dip_vectors(dip, dip_azimuth)
    u, v = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n), indexing="ij")
    pts = (
        np.asarray(center)
        + u.ravel()[:, None] * half * strike
        + v.ravel()[:, None] * half * down_dip
    )
    return pts, normal, down_dip


def add_fault(model, name, center, dip, dip_azimuth, displacement, *, points=None,
               profile=None, axes=(4.0, 4.0, 4.0), faults=None, steps=1, nelements=4000):
    plane, normal, down_dip = plane_points(center, dip, dip_azimuth)
    if points is None:
        points = plane
    inside = np.all((points > -0.1) & (points < 1.1), axis=1)
    points = points[inside]
    df = pd.DataFrame(
        {
            "X": points[:, 0], "Y": points[:, 1], "Z": points[:, 2], "val": 0.0,
            "nx": normal[0], "ny": normal[1], "nz": normal[2], "coord": 0,
            "feature_name": name,
        }
    )[COLUMNS]
    fault = model.create_and_add_fault(
        name,
        displacement,
        data=df,
        fault_slip_vector=-down_dip,  # restoration direction: normal fault
        fault_normal_vector=normal,
        fault_center=np.asarray(center, dtype=float),
        major_axis=axes[0],
        intermediate_axis=axes[1],
        minor_axis=axes[2],
        faultfunction=profile if profile is not None else constant_profile(),
        faults=faults or [],
        nelements=nelements,
        fault_buffer=0.5,
    )
    fault.steps = steps
    return fault


def add_layers(model, faults_oldest_first, n=800, seed=7, nelements=8000):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0.0, 1.0, size=(n, 3))
    restored = pts.copy()
    for f in reversed(faults_oldest_first):
        restored = f.apply_to_points(restored)
    df = pd.DataFrame(
        {
            "X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": restored[:, 2],
            "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": "strat",
        }
    )[COLUMNS]
    return model.create_and_add_foliation(
        "strat", data=df, faults=faults_oldest_first, nelements=nelements, buffer=0.3
    )


def interior_points(n=400, seed=3, lo=0.15, hi=0.85):
    return np.random.default_rng(seed).uniform(lo, hi, size=(n, 3))


YOUNG = dict(center=(0.65, 0.5, 0.5), dip=65, dip_azimuth=270, displacement=0.12)
OLD = dict(center=(0.3, 0.5, 0.4), dip=55, dip_azimuth=90, displacement=0.1)


def young_side(points):
    normal, _, _ = dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    return (np.asarray(points) - np.asarray(YOUNG["center"])) @ normal


def young_forward(points):
    """Move restored points forward through the young fault (exact: a translation)."""
    _, down_dip, _ = dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    points = np.array(points, dtype=float)
    hanging = young_side(points) > 0
    points[hanging] += YOUNG["displacement"] * down_dip
    return points
