from os.path import dirname, join
from pathlib import Path
import numpy as np
import pandas as pd


def load_claudius():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/claudius.csv')))
    bb = np.loadtxt(join(module_path, Path("data/claudiusbb.txt")))
    return data, bb


def load_noddy_single_fold():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/onefolddata.csv')))
    bb = np.loadtxt(join(module_path, Path('data/onefoldbb.txt')))
    return data, bb


def load_laurent2016():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/refolded_fold.csv')))
    bb = np.loadtxt(join(module_path, Path('data/refolded_bb.txt')))
    return data, bb

def load_duplex():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/duplex.csv')))
    bb = np.loadtxt(join(module_path, Path('data/duplexbb.txt')))
    return data, bb

def load_grose2017():
    pass


def load_grose2018():
    pass


def load_grose2019():
    pass


def load_intrusion():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path,Path('data/intrusion.csv')))
    bb = np.loadtxt(join(module_path,Path('data/intrusionbb.txt')))
    return data, bb
def load_unconformity():
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path,Path('data/unconformity.csv')))
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

