from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer

from LoopStructural.datasets import load_claudius #demo data

import pandas as pd
import glob
import numpy as np
import logging
logging.getLogger().setLevel(logging.INFO)

data, bb = load_claudius()#claudius.get_data()
bb[1,0]+=200
bb[0,0]-=200
bb[1,1]+=200
bb[0,1]-=200
bb[1,2]+=200
bb[0,2]-=200

model = GeologicalModel(bb[0,:],bb[1,:])
data['random'] = np.random.random(data.shape[0])
model.set_model_data(data[data['random'] < 0.01])#[np.isnan(data['val'])])
strati = model.create_and_add_foliation("strati",
                                           interpolatortype="surfe",
                                        )
