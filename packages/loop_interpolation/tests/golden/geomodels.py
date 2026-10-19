"""Synthetic geological models with an analytic scalar field.

Every model is a stratigraphic scalar field ``s(x, y, z)`` on the unit cube.
The horizons are the iso-surfaces ``s = level``. Data are sampled from the
analytic field, so the error of a solution can be measured against the truth.
All random numbers use fixed seeds: the data are identical on every run.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.optimize import brentq

HORIZONS = (0.3, 0.5, 0.7)


@dataclass
class Fold:
    """Cylindrical fold. Fold frame: gx = x (across), gy = y (along the axis, pointing to -y)."""

    axis: np.ndarray | None
    limb_angle_deg: Callable[[np.ndarray], np.ndarray]
    # Refolded folds: the axis is dgy rotated about -dgx by this angle (a function of gy).
    axis_angle_deg: Callable[[np.ndarray], np.ndarray] | None = None


@dataclass
class Case:
    name: str
    title: str
    why: str  # what a failure of this case means
    truth: Callable[[np.ndarray], np.ndarray]
    truth_grad: Callable[[np.ndarray], np.ndarray]
    values: np.ndarray  # (n, 4) x y z value
    normals: np.ndarray  # (n, 6) x y z nx ny nz
    tangents: np.ndarray = field(default_factory=lambda: np.zeros((0, 6)))
    fold: Fold | None = None
    engines: tuple[str, ...] = ("FD", "P1", "P2")
    # Corners of the region that accuracy is measured in.
    eval_box: tuple[tuple[float, float, float], tuple[float, float, float]] = (
        (0.1, 0.1, 0.1),
        (0.9, 0.9, 0.9),
    )
    # Axis normal to the figure section (1 = the XZ plane at y = 0.5).
    section_axis: int = 1


def _horizon_points(truth, xy: np.ndarray, levels=HORIZONS) -> np.ndarray:
    """Points on each horizon below the (x, y) locations (a borehole-style sample)."""
    out = []
    for x, y in xy:

        def f(z, level=None):
            return truth(np.array([[x, y, z]]))[0] - level

        for level in levels:
            if f(0.0, level) * f(1.0, level) < 0:
                z = brentq(f, 0.0, 1.0, args=(level,), xtol=1e-12)
                out.append([x, y, z, level])
    return np.array(out)


def _unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


# --------------------------------------------------------------------------- planar
def _plane(dip_deg=30.0, azimuth_deg=40.0):
    d, a = np.deg2rad(dip_deg), np.deg2rad(azimuth_deg)
    n = np.array([np.sin(d) * np.cos(a), np.sin(d) * np.sin(a), np.cos(d)])
    return (lambda p: p @ n), (lambda p: np.tile(n, (len(p), 1)))


def planar_dip() -> Case:
    s, g = _plane()
    rng = np.random.default_rng(1)
    values = _horizon_points(s, rng.uniform(0.1, 0.9, (6, 2)))
    npts = rng.uniform(0.2, 0.8, (3, 3))
    return Case(
        "planar_dip",
        "Planar dipping contacts (dip 30, 3 horizons, 6 boreholes)",
        "Linear truth. Every interpolator must reproduce it. A failure here is a basic bug "
        "(sign, scaling, weights, data placement).",
        s, g, values, np.c_[npts, g(npts)],
    )


def planar_extrapolation() -> Case:
    s, g = _plane()
    rng = np.random.default_rng(2)
    xy = np.c_[rng.uniform(0.05, 0.35, 5), rng.uniform(0.1, 0.9, 5)]
    values = _horizon_points(s, xy)
    npts = np.c_[rng.uniform(0.05, 0.35, 3), rng.uniform(0.1, 0.9, 3), rng.uniform(0.2, 0.8, 3)]
    return Case(
        "planar_extrapolation",
        "Planar contacts, data only in x < 0.35",
        "Tests regularisation as an extrapolator into a data-free region (x > 0.65). "
        "The regularisation should continue the plane, not flatten or curl it.",
        s, g, values, np.c_[npts, g(npts)],
    )


def planar_noisy() -> Case:
    case = planar_dip()
    rng = np.random.default_rng(3)
    values = case.values.copy()
    values[:, 2] += rng.normal(0, 0.01, len(values))  # 1% of the box in z
    normals = case.normals.copy()
    normals[:, 3:6] = normals[:, 3:6] + rng.normal(0, 0.05, (len(normals), 3))
    case.name = "planar_noisy"
    case.title = "Planar contacts with noisy data (z +-0.01, normals +-3 deg)"
    case.why = (
        "Checks smoothing: the solution should average the noise, not pass through every "
        "noisy point. Large error = regularisation too weak or data over-weighted."
    )
    case.values, case.normals = values, normals
    return case


# --------------------------------------------------------------------------- dome
def dome() -> Case:
    a, c, w = 0.2, 0.5, 0.2

    def s(p):
        r2 = (p[:, 0] - c) ** 2 + (p[:, 1] - c) ** 2
        return p[:, 2] - a * np.exp(-r2 / (2 * w**2))

    def g(p):
        r2 = (p[:, 0] - c) ** 2 + (p[:, 1] - c) ** 2
        e = a * np.exp(-r2 / (2 * w**2)) / w**2
        return np.c_[e * (p[:, 0] - c), e * (p[:, 1] - c), np.ones(len(p))]

    rng = np.random.default_rng(4)
    xy = np.vstack(
        [[0.5, 0.5], [0.35, 0.5], [0.65, 0.5], [0.5, 0.35], [0.5, 0.65], rng.uniform(0.1, 0.9, (4, 2))]
    )
    values = _horizon_points(s, xy)
    npts = np.array([[0.4, 0.4, 0.5], [0.6, 0.6, 0.5], [0.3, 0.7, 0.6], [0.8, 0.3, 0.4], [0.5, 0.5, 0.5]])
    return Case(
        "dome",
        "Gaussian dome (relief 0.2, sigma 0.2, 3 horizons, 9 boreholes)",
        "Smooth curved truth. Tests whether regularisation flattens the crest and whether "
        "gradient data bend the horizons correctly.",
        s, g, values, np.c_[npts, g(npts)],
    )


# --------------------------------------------------------------------------- folds
FOLD_AMP, FOLD_WL = 0.12, 0.6


def _fold_truth():
    k = 2 * np.pi / FOLD_WL

    def h(x):
        return FOLD_AMP * np.sin(k * x)

    def dh(x):
        return FOLD_AMP * k * np.cos(k * x)

    def s(p):
        return p[:, 2] - h(p[:, 0])

    def g(p):
        return np.c_[-dh(p[:, 0]), np.zeros(len(p)), np.ones(len(p))]

    def limb(gx):
        # The library rotates dgx by -theta about the fold axis (R transposed in
        # FoldEvent.get_deformed_orientation). With axis (0,-1,0), dgz = dgx x axis = (0,0,-1)
        # points against the data normals, which is the library convention (fold_norm = -1).
        # The fold direction must be ~ (1, 0, h'), which needs theta = -atan(h').
        return -np.rad2deg(np.arctan(dh(gx)))

    return s, g, Fold(axis=np.array([0.0, -1.0, 0.0]), limb_angle_deg=limb)


def fold_section() -> Case:
    """Dense cross-section data: plain interpolation of a cylindrical fold."""
    s, g, fold = _fold_truth()
    xy = np.c_[np.linspace(0.05, 0.95, 10), np.full(10, 0.5)]
    values = _horizon_points(s, xy)
    npts = np.c_[np.linspace(0.1, 0.9, 9), np.full(9, 0.5), np.full(9, 0.5)]
    return Case(
        "fold_section",
        "Cylindrical fold, dense section at y = 0.5 (no fold constraint)",
        "Plain interpolation of a fold from one cross-section. Along the axis (y) the "
        "regularisation must carry the section outward. Error far from y = 0.5 is expected; "
        "the golden file detects a change in the pattern.",
        s, g, values, np.c_[npts, g(npts)], fold=fold,
    )


def fold_sparse() -> Case:
    """Sparse data. Only the fold constraints can recover the fold shape."""
    s, g, fold = _fold_truth()
    values = _horizon_points(s, np.array([[0.2, 0.3], [0.8, 0.7]]), levels=(0.5,))
    p = np.array([[0.5, 0.5, 0.5]])
    return Case(
        "fold_sparse",
        "Cylindrical fold, 2 value points + 1 normal (fold constraints on)",
        "The Bermagui situation: little data, strong fold. The fold constraints must supply the "
        "gradient direction everywhere. The FD and P1 columns have no fold constraint and show "
        "what is gained. A drop in DFI/FDFold accuracy is a fold-constraint regression.",
        s, g, values, np.c_[p, g(p)], fold=fold, engines=("FD", "P1", "P2", "DFI", "FDFold"),
    )


def refolded_sparse() -> Case:
    """F2 folds (axis along y) refolded by F3 (axis along x): a dome-and-basin interference pattern.

    s = z - A sin(k2 x) - B sin(k3 y). The F2 hinge line plunges by atan(B k3 cos(k3 y)), so the
    fold axis rotates with the along-axis coordinate gy = y. Like Bermagui, the data are a few
    contacts and orientations, and the fold event supplies the rest.

    The fold event can only express the limb rotation as a function of gx, but the exact
    angle is tan(theta) = -h'(x) cos(plunge(y)). The event uses theta = -atan(h'), so the
    constraint direction is off by a few degrees where the plunge is steep. `model_floor_deg`
    in the tests measures this.
    """
    a2, wl2, b3, wl3 = 0.10, 0.5, 0.08, 0.8
    k2, k3 = 2 * np.pi / wl2, 2 * np.pi / wl3

    def s(p):
        return p[:, 2] - a2 * np.sin(k2 * p[:, 0]) - b3 * np.sin(k3 * p[:, 1])

    def g(p):
        return np.c_[-a2 * k2 * np.cos(k2 * p[:, 0]), -b3 * k3 * np.cos(k3 * p[:, 1]), np.ones(len(p))]

    def limb(gx):
        return -np.rad2deg(np.arctan(a2 * k2 * np.cos(k2 * gx)))

    def axis_rotation(gy):
        # +180: dgz = dgx x axis must point against the data normals (fold_norm = -1)
        return np.rad2deg(np.arctan(b3 * k3 * np.cos(k3 * gy))) + 180.0

    fold = Fold(axis=None, limb_angle_deg=limb, axis_angle_deg=axis_rotation)
    values = _horizon_points(s, np.array([[0.25, 0.25], [0.75, 0.5], [0.4, 0.8]]), levels=(0.5,))
    rng = np.random.default_rng(5)
    p = rng.uniform(0.15, 0.85, (6, 3))
    return Case(
        "refolded_sparse",
        "Refolded fold (F2 axis along y refolded by F3 about x): dome and basin, 3 contacts + 6 normals",
        "The Bermagui geometry in miniature: the hinge plunge changes along the axis, so a constant "
        "axis is wrong. Tests the fold_axis_rotation path of the fold interpolators (DFI, FDFold). "
        "FD, P1 and P2 have no fold constraint and show what the fold event adds.",
        s, g, values, np.c_[p, g(p)], fold=fold, engines=("FD", "P1", "P2", "DFI", "FDFold"),
    )


ALL = {f.__name__: f for f in (planar_dip, planar_extrapolation, planar_noisy, dome, fold_section, fold_sparse, refolded_sparse)}


def build(name: str) -> Case:
    return ALL[name]()
