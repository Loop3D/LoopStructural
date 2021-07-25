def test_import_model():
    success = True
    try:
        from LoopStructural import GeologicalModel
    except ImportError:
        success = False
    assert success == True

# def test_import_visualisation():
#     success = True    
#     try:
#         from LoopStructural.visualisation import LavaVuModelViewer
#     except ImportError:
#         success = False
#     assert success == True

# def test_import_rotation_angle_plotter():
#     success = True
#     try:
#         from LoopStructural.visualisation import RotationAnglePlotter
#     except ImportError:
#         success = False
#     assert success == True

def test_import_geological_feature():
    success = True
    try:
        from LoopStructural.modelling.features import GeologicalFeature
    except ImportError:
        success = False
    assert success == True