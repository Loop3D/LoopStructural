import numpy as np
import pandas as pd
import pytest
from loop_interpolation import DiscreteFoldInterpolator, P2FoldInterpolator

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold

data, boundary_points = load_noddy_single_fold()


def _model(interpolatortype=None, nelements=2000):
    mdata = pd.concat([data[:100], data[data["feature_name"] == "s1"]])
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    fold_frame = model.create_and_add_fold_frame("s1", nelements=nelements)
    kwargs = {} if interpolatortype is None else {"interpolatortype": interpolatortype}
    s0 = model.create_and_add_folded_foliation(
        "s0",
        fold_frame=fold_frame,
        nelements=nelements,
        fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        limb_wl=1,
        **kwargs,
    )
    model.update()
    return model, s0


def test_default_is_dfi():
    _, s0 = _model()
    assert isinstance(s0.builder.interpolator, DiscreteFoldInterpolator)


def test_unknown_type_falls_back_to_dfi():
    _, s0 = _model("FDI")
    assert isinstance(s0.builder.interpolator, DiscreteFoldInterpolator)


def test_p2fold_foliation_builds_and_evaluates():
    model, s0 = _model("P2Fold")
    interp = s0.builder.interpolator
    assert isinstance(interp, P2FoldInterpolator)
    assert interp.fold is s0.fold
    xyz = np.random.default_rng(0).uniform(
        boundary_points[0] + 1, boundary_points[1] - 1, size=(200, 3)
    )
    v = s0.evaluate_value(xyz)
    assert np.isfinite(v).all()
    assert np.isfinite(s0.evaluate_gradient(xyz)).all()


def test_p2fold_save_load(tmp_path):
    pytest.importorskip("dill")
    model, s0 = _model("P2Fold")
    xyz = np.random.default_rng(0).uniform(
        boundary_points[0] + 1, boundary_points[1] - 1, size=(200, 3)
    )
    v = s0.evaluate_value(xyz)
    path = tmp_path / "p2fold.pkl"
    model.to_file(str(path))
    loaded = GeologicalModel.from_file(str(path))
    assert np.allclose(loaded["s0"].evaluate_value(xyz), v, equal_nan=True)


def test_p2fold_folded_fold_frame():
    mdata = data
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(mdata)
    ff = model.create_and_add_fold_frame("s1", nelements=2000)
    frame = model.create_and_add_folded_fold_frame(
        "s0",
        fold_frame=ff,
        interpolatortype="P2Fold",
        nelements=2000,
        fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
        limb_wl=1,
    )
    assert isinstance(frame[0].builder.interpolator, P2FoldInterpolator)
