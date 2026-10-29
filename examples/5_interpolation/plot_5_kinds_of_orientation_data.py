"""
5e. Kinds of orientation data, with InterpolatorBuilder
=======================================================

Contact points tell the interpolator where a surface is. Orientation data tell it which way
the surface faces, but there are several ways to say that, and they carry different amounts of
information:

* a **tangent** is one direction that lies *in* the surface (for example a strike line or a
  lineation). One tangent leaves the surface free to rotate about it,
* **two tangents** that are not parallel fix the plane, and so the direction of the gradient,
* a **gradient** constraint gives the direction of the gradient (the surface normal) but not its
  length,
* a **normal** constraint gives the full gradient vector, direction and length. Its length
  says how quickly the scalar field changes across the surface, so it also fixes the scale of
  the field.

This example uses ``InterpolatorBuilder``, which assembles an interpolator from a bounding box
and arrays of constraints, and measures what each kind of orientation data adds when there are
few contact points. The truth is the same gaussian dome as in the earlier examples.
"""

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from LoopStructural import BoundingBox, InterpolatorBuilder

# %%
# The model and its data
# ----------------------

RELIEF, CENTRE, WIDTH = 0.2, 0.5, 0.2
HORIZONS = (0.3, 0.5, 0.7)


def truth(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    return xyz[:, 2] - RELIEF * np.exp(-r2 / (2 * WIDTH**2))


def truth_gradient(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    e = RELIEF * np.exp(-r2 / (2 * WIDTH**2)) / WIDTH**2
    return np.c_[e * (xyz[:, 0] - CENTRE), e * (xyz[:, 1] - CENTRE), np.ones(len(xyz))]


def borehole_contacts(xy):
    points = []
    for x, y in xy:
        for level in HORIZONS:

            def misfit(z, x=x, y=y, level=level):
                return truth(np.array([[x, y, z]]))[0] - level

            if misfit(0.0) * misfit(1.0) < 0:
                points.append([x, y, brentq(misfit, 0.0, 1.0), level])
    return np.array(points)


# %%
# Few contacts: two boreholes, so six contact points. Five orientation measurements, each
# written four ways.

contacts = borehole_contacts(np.array([[0.5, 0.5], [0.35, 0.5]]))
print(f"{len(contacts)} contact points")

orientation_xyz = np.array(
    [[0.4, 0.4, 0.5], [0.6, 0.6, 0.5], [0.3, 0.7, 0.6], [0.8, 0.3, 0.4], [0.5, 0.5, 0.5]]
)
gradient = truth_gradient(orientation_xyz)
unit_normal = gradient / np.linalg.norm(gradient, axis=1, keepdims=True)

# two directions in the surface at each point: ``first`` is horizontal where possible
reference = np.where(
    np.abs(unit_normal[:, [2]]) > 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 0, 1.0]])
)
first = np.cross(unit_normal, reference)
first /= np.linalg.norm(first, axis=1, keepdims=True)
second = np.cross(unit_normal, first)

constraints = {
    "contacts only": {},
    "one tangent per point": {"tangent": np.c_[orientation_xyz, first]},
    "two tangents per point": {
        "tangent": np.vstack([np.c_[orientation_xyz, first], np.c_[orientation_xyz, second]])
    },
    "gradient direction": {"gradient": np.c_[orientation_xyz, unit_normal]},
    "normal (direction and length)": {"normal": np.c_[orientation_xyz, gradient]},
}

# %%
# Building each interpolator
# --------------------------
# ``InterpolatorBuilder`` takes the interpolator type, the bounding box of the model and the
# approximate number of elements. Each ``add_*_constraints`` call returns the builder, so the
# calls chain, and ``build()`` returns the interpolator. The arrays are:
#
# ===============  ==========================================
# method           columns
# ===============  ==========================================
# value            x, y, z, value
# normal           x, y, z, nx, ny, nz
# gradient         x, y, z, gx, gy, gz (direction only)
# tangent          x, y, z, tx, ty, tz (lies in the surface)
# ===============  ==========================================

bounding_box = BoundingBox(origin=[0.0, 0.0, 0.0], maximum=[1.0, 1.0, 1.0])


def build(extra):
    builder = InterpolatorBuilder(
        "FDI", bounding_box=bounding_box, nelements=4000, buffer=0.1
    ).add_value_constraints(contacts)
    for kind, array in extra.items():
        getattr(builder, f"add_{kind}_constraints")(array)
    return builder.build()


interpolators = {label: build(extra) for label, extra in constraints.items()}

# %%
# What each kind of data adds
# ---------------------------
# Compare with the truth on a regular set of interior points.

axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T


def angle_degrees(a, b):
    cosine = np.sum(a * b, axis=1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))


rows = {}
for label, interpolator in interpolators.items():
    error = interpolator.evaluate_value(points) - truth(points)
    angle = angle_degrees(interpolator.evaluate_gradient(points), truth_gradient(points))
    rows[label] = {
        "rms error": np.sqrt(np.nanmean(error**2)),
        "median gradient angle": np.nanmedian(angle),
    }
print(pd.DataFrame(rows).T.round(3).to_string())

# %%
# What to take from this
# ----------------------
# * Each step adds information and the error falls in order: contacts only, one tangent, two
#   tangents or a gradient direction, then a full normal.
# * One tangent helps only a little. It removes one of the two directions the surface could
#   rotate in, and leaves the other free.
# * Two non-parallel tangents are nearly equivalent to a gradient direction: both fix the
#   plane, and the errors are within a few percent of each other.
# * A normal is the most informative because its length constrains the scale of the field as
#   well as its direction. That is also why it is the usual choice in ``GeologicalModel``
#   (the ``nx, ny, nz`` columns); use a gradient direction when you have no reliable measure of
#   how fast the field changes, and tangents when the data really are lines in the surface.
#
# These are errors from one set of six contacts and five orientations. With more contact points
# the differences shrink, because the contacts alone already constrain the surface well.
