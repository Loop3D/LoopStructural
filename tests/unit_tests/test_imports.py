def test_all_imports():
    success = True
    try:
        from LoopStructural import GeologicalModel
        from LoopStructural.visualisation import LavaVuModelViewer
        from LoopStructural.visualisation import RotationAnglePlotter
        from LoopStructural.utils import process_map2loop, build_model
        from LoopStructural.modelling.features import GeologicalFeature
    except ImportError:
        success = False
    assert success == True