import warnings


def test_datatypes_reexports_geometry_symbols():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from LoopStructural.datatypes import (
            BoundingBox,
            Surface,
            ValuePoints,
            VectorPoints,
        )

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)

    from LoopStructural import geometry

    assert BoundingBox is geometry.BoundingBox
    assert Surface is geometry.Surface
    assert ValuePoints is geometry.ValuePoints
    assert VectorPoints is geometry.VectorPoints
