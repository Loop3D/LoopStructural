from LoopStructural import GeologicalModel
from LoopStructural.modelling.fault import FaultSegment
import pandas as pd
def test_create_and_add_fault():
    model = GeologicalModel([0,0,0],[1,1,1])
    data = pd.DataFrame([[.5,.5,.5,0,1,0,0,'fault',0],
                         [.5,.5,.5,1,0,0,1,'fault',0],
                         [.5,.5,.5,0,0,1,2,'fault',0],],
                         columns=['X','Y','Z','nx','ny','nz','coord','feature_name','val'])
    model.data = data
    model.create_and_add_fault('fault',1,nelements=1e4,force_mesh_geometry=True)
    assert isinstance(model['fault'],FaultSegment)