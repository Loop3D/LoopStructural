"""
Mt Bruce
========
This example shows how to use the output from map2loop to create a 3D geological model
of the Mt Bruce area in the Hamersley. The dataset has been compiled into a pickle object/
"""
from LoopStructural import GeologicalModel
from LoopStructural.modelling.fault.fault_function import BaseFault
from LoopStructural.visualisation import LavaVuModelViewer
import pandas as pd
import numpy as np
import pickle


# m2l_data = pickle.load(open('mt_bruce.pkl','rb'))