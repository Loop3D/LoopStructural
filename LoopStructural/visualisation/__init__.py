try:
    from loopstructuralvisualisation import (
        Loop3DView,
        RotationAnglePlotter,
        Loop2DView,
        StratigraphicColumnView,
    )
except ImportError as e:
    print("Please install the loopstructuralvisualisation package")
    print("pip install loopstructuralvisualisation")
    raise e
