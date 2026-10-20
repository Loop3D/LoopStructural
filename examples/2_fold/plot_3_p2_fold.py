"""
2c. Quadratic fold interpolation
==============================

``P2Fold`` uses quadratic tetrahedra and samples the fold constraints inside
each tetrahedron. Select it when constructing a folded foliation or the first
coordinate of a folded fold frame. ``gauss4`` is the default sample rule.

Start with a small mesh and inspect the realised unknown count. The current
factory interprets ``nelements`` as a grid-cell budget, then splits each cell
into five tetrahedra; P2 also adds nodes on the tetrahedron edges.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold

# This bundled synthetic dataset makes the example independent of Bermagui data.
data, bounds = load_noddy_single_fold()
model = GeologicalModel(bounds[0], bounds[1])
model.set_model_data(pd.concat([data[:100], data[data["feature_name"] == "s1"]]))
frame = model.create_and_add_fold_frame("s1", nelements=1000, buffer=0.2)
folded = model.create_and_add_folded_foliation(
    "s0",
    fold_frame=frame,
    interpolatortype="P2Fold",
    nelements=1000,
    buffer=0.0,
    fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
    limb_wl=1,
    fold_weights={"sample_rule": "gauss4"},
    solver="cg",
    solver_kwargs={"maxiter": 20000, "rtol": 1e-8},
)

# Builds are lazy: evaluation triggers the fit and solve.
x = np.linspace(
    bounds[0, 0] + 0.1 * np.ptp(bounds[:, 0]), bounds[1, 0] - 0.1 * np.ptp(bounds[:, 0]), 60
)
z = np.linspace(
    bounds[0, 2] + 0.1 * np.ptp(bounds[:, 2]), bounds[1, 2] - 0.1 * np.ptp(bounds[:, 2]), 60
)
xx, zz = np.meshgrid(x, z)
points = np.c_[xx.ravel(), np.full(xx.size, np.mean(bounds[:, 1])), zz.ravel()]
values = folded.evaluate_value(points)
assert np.isfinite(values).all()
interpolator = folded.builder.interpolator
print(f"P2Fold: {interpolator.support.n_elements:,} tetrahedra, {interpolator.dof:,} unknowns")

fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
contours = ax.contourf(xx, zz, values.reshape(xx.shape), levels=20, cmap="viridis")
ax.contour(xx, zz, values.reshape(xx.shape), levels=10, colors="white", linewidths=0.5)
ax.set(xlabel="X", ylabel="Z", title="P2Fold with gauss4: section through the synthetic fold")
fig.colorbar(contours, ax=ax, label="Folded scalar field (arbitrary units)")
plt.show()
