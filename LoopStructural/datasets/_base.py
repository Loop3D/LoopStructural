from os.path import dirname, join
from pathlib import Path
import numpy as np
import pandas as pd


def load_claudius():
    """Model dataset sampled from 3D seismic data


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/claudius.csv')))
    bb = np.loadtxt(join(module_path, Path("data/claudiusbb.txt")))
    return data, bb


def load_noddy_single_fold():
    """Model dataset for plunging cylindrical fold 


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """

    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/onefolddata.csv')))
    bb = np.loadtxt(join(module_path, Path('data/onefoldbb.txt')))
    return data, bb


def load_laurent2016():
    """Model dataset for refolded fold 


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/refolded_fold.csv')))
    bb = np.loadtxt(join(module_path, Path('data/refolded_bb.txt')))
    return data, bb

def load_duplex():
    """Model dataset for synthetic duplex example


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path, Path('data/duplex.csv')))
    bb = np.loadtxt(join(module_path, Path('data/duplexbb.txt')))
    return data, bb

def load_grose2017():
    """Model dataset for Cape Conran


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    pass


def load_grose2018():
    """Model dataset for synthetic parasitic fold series


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    pass


def load_grose2019():
    """Model dataset for Davenport ranges


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    pass


def load_intrusion():
    """Model dataset for a faulted intrusion


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
    module_path = dirname(__file__)
    data = pd.read_csv(join(module_path,Path('data/intrusion.csv')))
    bb = np.loadtxt(join(module_path,Path('data/intrusionbb.txt')))
    return data, bb
def load_unconformity():
    """Model dataset sampled for a model containing an unconformity


    Returns
    -------
    tuple
        pandas data frame with loopstructural dataset and numpy array for bounding box
    """
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

