import importlib.util


def test_import_model():

    success = True
    if importlib.util.find_spec("LoopStructural", 'GeologicalModel') is None:
        success = False
    assert success


# def test_import_visualisation():
#     success = True
#     try:
#         from LoopStructural.visualisation import Loop3DView
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


# def test_import_geological_feature():

#     success = True
#     if importlib.util.find_spec('LoopStructural.modelling.features', ''):
#         success = False
#     assert success
