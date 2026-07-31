"""
1i. Performance of a deep unconformity stack
=============================================
This example builds a stack of 10 boundaries alternating between
**erosional** unconformities (:code:`add_unconformity`) and **onlap**
unconformities (:code:`add_onlap_unconformity`), with a fault inserted
partway up the stack, and times how long it takes to evaluate the scalar
field of the *oldest* feature in the stack.

Each boundary that is added should only affect features added *after* it -
older features should never need to know about younger boundaries. This
example also reproduces the pre-fix behaviour of
:code:`add_onlap_unconformity`, where the backward search for existing
features to attach the onlap region to did not stop at the previous
unconformity, so *every* older feature (all the way back to the oldest
one in the model) ended up carrying regions from onlap surfaces added much
later - including one on the far side of the fault. That made evaluating
the oldest feature dramatically slower than it needed to be, because
evaluating those spurious regions also evaluated the fault restoration.
"""

import time
import types

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import FeatureType, UnconformityFeature

######################################################################
# Reproducing the pre-fix behaviour
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ``_buggy_add_onlap_unconformity`` is a copy of ``add_onlap_unconformity``
# as it existed before the fix: it ``continue``\ s past existing
# unconformities instead of ``break``\ ing, so it keeps walking all the way
# back through the entire feature history rather than stopping at the last
# boundary.


def _buggy_add_onlap_unconformity(self, feature, value, index=None):
    feature.regions = []
    uc_feature = UnconformityFeature(feature, value, False, onlap=True)
    feature.add_region(uc_feature.inverse())
    for f in reversed(self.features):
        if f.type == FeatureType.UNCONFORMITY:
            continue
        if f.type == FeatureType.FAULT:
            continue
        if f != feature:
            f.add_region(uc_feature)
    self._add_feature(uc_feature.inverse(), index=index)
    return uc_feature


######################################################################
# Building the stack
# ~~~~~~~~~~~~~~~~~~~
# 11 units are separated by 10 boundaries (alternating erosional/onlap), with
# a fault inserted early in the sequence - well before most of the onlap
# boundaries are added. ``unit0`` is the oldest feature in the model.

N_BOUNDARIES = 10


def build_stacked_model(buggy_onlap: bool) -> GeologicalModel:
    unit_names = [f"unit{i}" for i in range(N_BOUNDARIES + 1)]

    rows = [[100, 100, 20 + i * 15, 0, 0, 1, 0, name] for i, name in enumerate(unit_names)]
    rows.append([700, 100, 190, 1, 0, 0, np.nan, "fault"])
    data = pd.DataFrame(rows, columns=["X", "Y", "Z", "nx", "ny", "nz", "val", "feature_name"])

    model = GeologicalModel(np.zeros(3), np.array([1000, 1000, 200]))
    model.data = data

    if buggy_onlap:
        model.add_onlap_unconformity = types.MethodType(_buggy_add_onlap_unconformity, model)

    model.create_and_add_foliation(unit_names[0], buffer=0.0)
    for i in range(N_BOUNDARIES):
        if i % 2 == 0:
            model.add_unconformity(model[unit_names[i]], 0)
        else:
            model.add_onlap_unconformity(model[unit_names[i]], 0)
        if i == 1:
            # insert a fault early in the stack - only features added from
            # here onwards should ever need to restore points through it
            model.create_and_add_fault(
                "fault",
                50,
                minor_axis=300,
                major_axis=500,
                intermediate_axis=300,
                fault_center=[700, 500, 0],
            )
        model.create_and_add_foliation(unit_names[i + 1], buffer=0.0)

    return model, unit_names[0]


######################################################################
# Timing the oldest feature
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# The same 100x100 grid used elsewhere in the unconformity/fault examples is
# evaluated repeatedly for the oldest feature ("unit0") in both the fixed and
# the (reproduced) buggy model.

xx, zz = np.meshgrid(np.linspace(0, 1000, 100), np.linspace(0, 200, 100))
yy = np.zeros_like(xx) + 500
points = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T


def time_evaluation(model, feature_name, n=10):
    feature = model[feature_name]
    feature.evaluate_value(points)  # warm-up / build
    t0 = time.perf_counter()
    for _ in range(n):
        feature.evaluate_value(points)
    t1 = time.perf_counter()
    return (t1 - t0) / n * 1e3


fixed_model, oldest_name = build_stacked_model(buggy_onlap=False)
buggy_model, _ = build_stacked_model(buggy_onlap=True)

fixed_ms = time_evaluation(fixed_model, oldest_name)
buggy_ms = time_evaluation(buggy_model, oldest_name)

print(f"Stack of {N_BOUNDARIES} alternating erosional/onlap unconformities + 1 fault")
print(f"Evaluating the oldest feature ('{oldest_name}') at {points.shape[0]} points:")
print(f"  fixed add_onlap_unconformity : {fixed_ms:8.3f} ms")
print(f"  buggy add_onlap_unconformity : {buggy_ms:8.3f} ms")
print(f"  speedup                      : {buggy_ms / fixed_ms:6.1f}x")

fig, ax = plt.subplots(figsize=(4, 4))
ax.bar(["fixed", "buggy (pre-fix)"], [fixed_ms, buggy_ms], color=["tab:green", "tab:red"])
ax.set_ylabel("mean evaluate_value time (ms)")
ax.set_title(f"Evaluating oldest feature '{oldest_name}'\nin a {N_BOUNDARIES}-boundary stack")
plt.tight_layout()
plt.show()
