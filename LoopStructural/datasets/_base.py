from os.path import dirname, join
import pandas as pd
import numpy as np

def load_claudius():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path,'data/claudius.pkl'))
    bb = np.loadtxt(join(module_path,"data/claudiusbb.txt"))
    return data, bb


def load_noddy_single_fold():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path,'data/onefolddata.pkl'))
    bb = np.loadtxt(join(module_path,'data/onefoldbb.txt'))
    return data, bb
def load_laurent2016():
    pass

def load_grose2017():
    pass

def load_grose2018():
    pass

def load_grose2019():
    pass

def load_intrusion():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path,'data/intrusion.pkl'))
    bb = np.loadtxt(join(module_path,'data/intrusiondbb.txt'))
    return data, bb