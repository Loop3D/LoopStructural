from os.path import dirname, join

import numpy as np
import pandas as pd


def load_claudius():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path, 'data/claudius.pkl'))
    bb = np.loadtxt(join(module_path, "data/claudiusbb.txt"))
    return data, bb


def load_noddy_single_fold():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path, 'data/onefolddata.pkl'))
    bb = np.loadtxt(join(module_path, 'data/onefoldbb.txt'))
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
    bb = np.loadtxt(join(module_path,'data/intrusionbb.txt'))
    return data, bb
def load_unconformity():
    module_path = dirname(__file__)
    data = pd.read_pickle(join(module_path,'data/unconformity.pkl'))
    bb = np.array([[0,0,0],
                   [4,6,4]]
                  )
    return data,bb
def value_headers():
    return ['X','Y','Z','val']

def strike_dip_headers():
    return ['X','Y','Z','strike','dip']

def normal_vector_headers():
    return ['X','Y','Z','nx','ny','nz']

