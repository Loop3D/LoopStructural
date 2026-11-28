"""Public import smoke tests; plotting and notebook behavior have separate suites."""

import loopstructuralvisualisation


def test_import():
    assert loopstructuralvisualisation is not None


def test_public_api_importable():
    from loopstructuralvisualisation import (
        Loop2DView,
        Loop3DView,
        RotationAnglePlotter,
        StratigraphicColumnView,
    )

    assert Loop2DView is not None
    assert Loop3DView is not None
    assert RotationAnglePlotter is not None
    assert StratigraphicColumnView is not None


def test_api_module_importable():
    from loopstructuralvisualisation.api import plot_block_model, plot_surface

    assert callable(plot_block_model)
    assert callable(plot_surface)
