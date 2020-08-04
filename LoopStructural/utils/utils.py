import logging

import numpy as np
import re

logger = logging.getLogger(__name__)


def strike_symbol(strike):
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))

    vec = np.array([0, 1])
    rotated = R @ vec
    vec2 = np.array([-0.5, 0])
    r2 = R @ vec2
    return rotated, r2
def get_levels():
    """dict for converting to logger levels from string


    Returns
    -------
    dict
        contains all strings with corresponding logging levels.
    """
    return {'info':logging.INFO,'warning':logging.WARNING,'error':logging.ERROR,'debug':logging.DEBUG}

def log_to_file(filename,level='info'):
    """Set the logging parameters for log file


    Parameters
    ----------
    filename : string
        name of file or path to file
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level,logging.WARNING)
    logging.basicConfig(level=level,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=filename,
                    filemode='w')

def log_to_console(level='warning'):
    """Set the level of logging to the console


    Parameters
    ----------
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level,logging.WARNING)

    changed_level = False
    for h in logging.getLogger().handlers:
        if type(h) is logging.StreamHandler:
            h.setLevel(level)
            changed_level = True
    if not changed_level:
        console = logging.StreamHandler()
        console.setLevel(level)
        # add the handler to the root logger
        logging.getLogger().addHandler(console)

def read_voxet(voxetname,propertyfile):
    """
    Read a gocad property file and the geometry information from the .vo file
    voxetname - is the path to the voxet file
    propertyfile is the path to the binary file
    Returns
    origin numpy array
    voxet_extent - is the length of each axis of the voxet
    N is the number of steps in the voxet
    array is the property values
    steps is the size of the step vector for the voxet 
    """
    array = np.fromfile(propertyfile,dtype='float32')
    array = array.astype('<f4') # little endian
    with open(voxetname,'r') as file:
        for l in file:
            if 'AXIS_O ' in l:
                origin = np.array(re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+",l)).astype(float)
            if 'AXIS_U ' in l:
                U = float(re.findall(r'[\d\.\d]+',l)[0])
            if 'AXIS_V ' in l:
                V = float(re.findall(r'[\d\.\d]+',l)[1])
            if 'AXIS_W ' in l:
                W = float(re.findall(r'[\d\.\d]+',l)[2])
            if 'AXIS_N ' in l:
                N = np.array(re.findall(r'[\d\.\d]+',l)).astype(int) 
    voxet_extent = np.array([U,V,W])
    steps = (voxet_extent ) / (N-1)
    return origin, voxet_extent, N, array, steps

def write_property_to_gocad_voxet(propertyfilename, propertyvalues):
    """
    This function writes a numpy array into the right format for a gocad
    voxet property file. This assumet there is a property already added to the .vo file,
    and is just updating the file.
    propertyfile - string giving the path to the file to write
    propertyvalues - numpy array nz,ny,nx ordering and in float format
    """
    propertyvalues = propertyvalues.astype('>f4') #big endian
#     array = propertyvalues.newbyteorder()
    propertyvalues.tofile(propertyfilename)
    