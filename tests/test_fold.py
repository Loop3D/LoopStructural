from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold

import pandas as pd

data, boundary_points = load_noddy_single_fold()
data.head()
def test_average_fold_axis():
    mdata = pd.concat([data[:100],data[data['feature_name']=='s1']])
    model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame('s1',nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation('s0',
                                                   fold_frame,
                                                    nelements=10000,
                                                    av_fold_axis=True
                                                   # fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                                   # limb_wl=1
                                                         )
def test_fixed_fold_axis():
    mdata = pd.concat([data[:100],data[data['feature_name']=='s1']])
    model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame('s1',nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation('s0',
                                                   fold_frame,
                                                    nelements=10000,
                                                    # av_fold_axis=True
                                                   fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                                   # limb_wl=1
                                                         )

def test_fixed_wavelength():
    mdata = pd.concat([data[:100],data[data['feature_name']=='s1']])
    model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame('s1',nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation('s0',
                                                   fold_frame,
                                                    nelements=10000,
                                                    # av_fold_axis=True
                                                   fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                                   limb_wl=1
                                                         )

def test_no_fold_frame():
    mdata = pd.concat([data[:100],data[data['feature_name']=='s1']])
    model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame('s1',nelements=10000)
    stratigraphy = model.create_and_add_folded_foliation('s0',
                                                   # fold_frame,
                                                    nelements=10000,
                                                    # av_fold_axis=True
                                                   fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
                                                   limb_wl=1
                                                         )