try:
    from loopstructuralvisualisation import (
        Loop2DView,
        Loop3DView,
        RotationAnglePlotter,
        StratigraphicColumnView,
    )
except ImportError as e:
    print("Please install the loopstructuralvisualisation package")
    print("pip install loopstructuralvisualisation")
    raise e
