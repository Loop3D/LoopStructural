import numpy as np
import re
from ..utils import getLogger

logger = getLogger(__name__)


def strike_symbol(strike):
    """Create rotation vectors for geological strike symbols.

    Generate rotation matrix and vectors for displaying geological strike symbols
    based on the strike angle.

    Parameters
    ----------
    strike : float
        Strike angle in degrees

    Returns
    -------
    rotated : np.ndarray
        Rotated vector for primary strike direction
    r2 : np.ndarray
        Rotated vector for secondary strike direction
    """
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


def read_voxet(voxetname, propertyfile):
    """Read a GOCAD property file and the geometry information from the .vo file.

    Parameters
    ----------
    voxetname : str
        Path to the voxet (.vo) file
    propertyfile : str
        Path to the binary property file

    Returns
    -------
    origin : np.ndarray
        Origin point of the voxet as numpy array
    voxet_extent : np.ndarray
        Length of each axis of the voxet
    N : np.ndarray
        Number of steps in each direction of the voxet
    array : np.ndarray
        Property values from the binary file
    steps : np.ndarray
        Size of the step vector for the voxet
    """
    array = np.fromfile(propertyfile, dtype="float32")
    array = array.astype("<f4")  # little endian
    with open(voxetname, "r") as file:
        for l in file:
            if "AXIS_O " in l:
                origin = np.array(re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", l)).astype(float)
            if "AXIS_U " in l:
                U = float(re.findall(r"[\d\.\d]+", l)[0])
            if "AXIS_V " in l:
                V = float(re.findall(r"[\d\.\d]+", l)[1])
            if "AXIS_W " in l:
                W = float(re.findall(r"[\d\.\d]+", l)[2])
            if "AXIS_N " in l:
                N = np.array(re.findall(r"[\d\.\d]+", l)).astype(int)
    voxet_extent = np.array([U, V, W])
    steps = (voxet_extent) / (N - 1)
    return origin, voxet_extent, N, array, steps


def write_property_to_gocad_voxet(propertyfilename, propertyvalues):
    """Write a numpy array to a GOCAD voxet property file format.

    This function writes a numpy array into the right format for a GOCAD
    voxet property file. This assumes there is a property already added to 
    the .vo file, and is just updating the file.

    Parameters
    ----------
    propertyfilename : str
        Path to the file where the property values will be written
    propertyvalues : np.ndarray
        Numpy array with nz,ny,nx ordering and in float format

    Notes
    -----
    The property values are converted to big-endian format before writing.
    """
    propertyvalues = propertyvalues.astype(">f4")  # big endian
    #     array = propertyvalues.newbyteorder()
    propertyvalues.tofile(propertyfilename)
