"""The ``solver`` build argument must reach the interpolator before setup."""

import logging

import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius


def _claudius_model(**build_kwargs):
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    feature = model.create_and_add_foliation(
        "strati", interpolatortype="FDI", nelements=2000, **build_kwargs
    )
    return model, feature


def test_build_sets_the_interpolator_solver():
    model, feature = _claudius_model(solver="lsmr")
    model.update(progressbar=False)
    assert feature.builder.interpolator.solver == "lsmr"


def test_native_solver_set_in_build_arguments_builds_no_rows(caplog):
    pytest.importorskip("loop_interpolation_core")
    model, feature = _claudius_model(solver="native")
    with caplog.at_level(logging.WARNING):
        model.update(progressbar=False)
    interpolator = feature.builder.interpolator
    assert interpolator.solver == "native"
    assert interpolator.latest_native_report["converged"]
    # A native setup stores the data but builds no sparse rows.
    assert len(interpolator.constraints) == 0
    assert "saves memory only when" not in caplog.text
