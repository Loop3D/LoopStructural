from LoopStructural.modelling.features import StructuralFrame, StructuralFrameBuilder, GeologicalFeature
from LoopStructural import GeologicalModel
def test_structural_frame():
    coordinate_0 = GeologicalFeature('coord0',None)
    coordinate_1 = GeologicalFeature('coord1',None)
    coordinate_2 = GeologicalFeature('coord2',None)
    frame = StructuralFrame('structural_frame',[coordinate_0,coordinate_1,coordinate_2])
    assert (frame != None)
    assert frame.name == 'structural_frame'
def set_model():
    model = GeologicalModel(np.zeros(),np.ones())
    coordinate_0 = GeologicalFeature('coord0',None)
    coordinate_1 = GeologicalFeature('coord1',None)
    coordinate_2 = GeologicalFeature('coord2',None)
    frame = StructuralFrame('structural_frame',[coordinate_0,coordinate_1,coordinate_2])
    frame.set_model(model)
    assert frame.model == model
    assert frame[0].model == model
    assert frame[1].model == model
    assert frame[2].model == model

def get_item():
    coordinate_0 = GeologicalFeature('coord0',None)
    coordinate_1 = GeologicalFeature('coord1',None)
    coordinate_2 = GeologicalFeature('coord2',None)
    frame = StructuralFrame('structural_frame',[coordinate_0,coordinate_1,coordinate_2])
    assert frame[0] == coordinate_0
    assert frame[1] == coordinate_0
    assert frame[2] == coordinate_0

