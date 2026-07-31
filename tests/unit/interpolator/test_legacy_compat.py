import importlib

import pytest
from loop_interpolation._discrete_interpolator import (
    DiscreteInterpolator as LoopDiscreteInterpolator,
)
from loop_interpolation._geological_interpolator import (
    GeologicalInterpolator as LoopGeologicalInterpolator,
)
from loop_interpolation._operator import Operator as LoopOperator
from loop_interpolation._p1interpolator import P1Interpolator as LoopP1Interpolator

from LoopStructural.interpolators import (
    DiscreteInterpolator,
    GeologicalInterpolator,
    Operator,
    P1Interpolator,
)


def test_public_api_reexports_loop_interpolation_classes():
    assert issubclass(DiscreteInterpolator, LoopDiscreteInterpolator)
    assert issubclass(GeologicalInterpolator, LoopGeologicalInterpolator)
    assert issubclass(P1Interpolator, LoopP1Interpolator)
    assert issubclass(Operator, LoopOperator)


def test_legacy_module_paths_are_removed():
    for module_name in (
        "LoopStructural.interpolators._discrete_interpolator",
        "LoopStructural.interpolators._geological_interpolator",
        "LoopStructural.interpolators._operator",
        "LoopStructural.interpolators._p1interpolator",
    ):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(module_name)
