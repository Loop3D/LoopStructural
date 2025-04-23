import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.modelling.features import GeologicalFeature

data, boundary_points = load_noddy_single_fold()
data.head()


def test_average_fold_axis():
    mdata = pd.concat([data[:100], data[data["feature_name"] == "s1"]])
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame("s1", nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation(
        "s0",
        fold_frame,
        nelements=10000,
        av_fold_axis=True,
        # fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        # limb_wl=1
    )
    model.update()

    assert np.isclose(
        stratigraphy.fold.fold_axis,
        np.array([-6.51626577e-06, -5.00013645e-01, -8.66017526e-01]),
        rtol=1e-3,
        atol=1e-3,
    ).all()


def test_fixed_fold_axis():
    mdata = pd.concat([data[:100], data[data["feature_name"] == "s1"]])
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame("s1", nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation(
        "s0",
        fold_frame,
        nelements=10000,
        # av_fold_axis=True
        fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        # limb_wl=1
    )
    model.update()
    assert np.isclose(
        stratigraphy.fold.fold_axis, np.array([-6.51626577e-06, -5.00013645e-01, -8.66017526e-01])
    ).all()


def test_fixed_wavelength():
    mdata = pd.concat([data[:100], data[data["feature_name"] == "s1"]])
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame("s1", nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation(
        "s0",
        fold_frame,
        nelements=10000,
        # av_fold_axis=True
        fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        limb_wl=1,
    )
    model.update()
    assert isinstance(stratigraphy, GeologicalFeature)


def test_no_fold_frame():
    mdata = pd.concat([data[:100], data[data["feature_name"] == "s1"]])
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame("s1", nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation(
        "s0",
        # fold_frame,
        nelements=10000,
        # av_fold_axis=True
        fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        limb_wl=1,
    )
    model.update()
    assert isinstance(stratigraphy, GeologicalFeature)
    assert stratigraphy.fold.foldframe.name == "s1"
    assert fold_frame.name == "s1"
