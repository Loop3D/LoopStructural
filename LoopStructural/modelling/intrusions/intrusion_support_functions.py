## Support Functions for intrusion network simulated as the shortest path, and for simulations in general
import numpy as np
import pandas as pd
import random
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


def array_from_coords(df, section_axis, df_axis):
    """
    Create numpy array representing a section of the model
    from a dataframe containing coordinates and values

    Parameters
    ----------
    df: pandas dataframe, should have at least ['X', 'Y', 'Z', 'val_1'] columns
    section_axis: string 'X' or 'Y', the cross section represented by the arrays is along the section_axis
    df_axis: number of the column where the value on interest is

    Returns
    -------
    array: numpy array, contains values (e.g. velocities at each point, scalar field value at each point, etc)

    """

    if section_axis == "X":
        other_axis = "Y"
    elif section_axis == "Y":
        other_axis = "X"

    col = len(df.columns)
    if col < df_axis:
        logger.error("Finding shortest path, dataframe axis out of range")

    else:
        df.sort_values([other_axis, "Z"], ascending=[True, False], inplace=True)
        xys = df[other_axis].unique()
        zs = df["Z"].unique()
        rows = len(zs)
        columns = len(xys)
        array = np.zeros([rows, columns])
        n = 0
        for j in range(columns):
            for i in range(rows):
                array[i, j] = df.iloc[i + n, df_axis]
            n = n + rows
    return array


def find_inout_points(velocity_field_array, velocity_parameters):
    """
    Looks for the indexes of the inlet and outle in an array.
    Velocity parameters of anisotropies are used to find the indexes of inlet and outlet.
    It is assumed that velocity_parameter[0] correspond to the inlet anisotropy and
    velocity_parameter[len(velocity_parameter)-1] corresponds to the outlet anisotropy

    Parameters
    ----------
    velocity_field_array: numpy array, containing values of the velocity field used to find the shortest path
    velocity_parameters: list of numbers, each value correspond to a velocity assign to an anisotropy involved in intrusion emplacement

    Returns
    -------
    inlet: list of indexes, [row index in array, column index in array]
    outlet: list of indexes, [row index in array, column index in array]

    """
    inlet_point = [0, 0]
    outlet_point = [0, 0]

    inlet_velocity = velocity_parameters[0] + 0.1
    outlet_velocity = velocity_parameters[len(velocity_parameters) - 1] + 0.1

    k = 0
    for i in range(len(velocity_field_array[0])):
        if k == 1:
            break

        where_inlet_i = np.where(velocity_field_array[:, i] == inlet_velocity)

        if len(where_inlet_i[0]) > 0:
            inlet_point[0] = where_inlet_i[0][len(where_inlet_i[0]) - 1]
            inlet_point[1] = i
            k = 1
        else:
            continue

    k = 0
    for i in range(len(velocity_field_array[0])):
        i_ = len(velocity_field_array[0]) - 1 - i
        if k == 1:
            break

        where_outlet_i = np.where(velocity_field_array[:, i_] == outlet_velocity)
        if len(where_outlet_i[0]) > 0:
            outlet_point[0] = where_outlet_i[0][0]
            outlet_point[1] = i_
            k = 1
        else:
            continue

    return inlet_point, outlet_point


def grid_from_array(array, fixed_coord, lower_extent, upper_extent):

    """
    Create an numpy matrix of [i,j,x,y,z,values in array], given an array of 2 dimensions (any combination between x, y an z)

    Parameters
    ----------
    array: numpy array, two dimension. Represents a cross section of the model, and its values could be any property
    fixed_coord: list, containing coordinate and value,
            ie, [0,2] means section is in x=2, or [1, .45] means sections is in y= 0.45
            the cross section is along this coordinate
    lower_extent: numpy array 1x3, lower extent of the model
    upper_extent: numpy array 1x3, upper extent of the model

    Returns
    -------
    values: numpy matrix of [i,j,x,y,z,values in array]
            (i,j) indexed in array
            (x,y,z) coordinates considering lower and upper extent of model
            values, from array

    """

    spacing_i = len(array)  # number of rows
    spacing_j = len(array[0])  # number of columns
    values = np.zeros([spacing_i * spacing_j, 6])
    if fixed_coord[0] == "X":
        y = np.linspace(lower_extent[1], upper_extent[1], spacing_j)
        z = np.linspace(lower_extent[2], upper_extent[2], spacing_i)
        l = 0
        for j in range(spacing_j):
            for i in range(spacing_i):
                values[l] = [
                    i,
                    j,
                    fixed_coord[1],
                    y[j],
                    z[i],
                    array[spacing_i - 1 - i, j],
                ]
                l = l + 1

    if fixed_coord[0] == "Y":
        x = np.linspace(lower_extent[0], upper_extent[0], spacing_j)
        z = np.linspace(lower_extent[2], upper_extent[2], spacing_i)
        l = 0
        for j in range(spacing_j):
            for i in range(spacing_i):
                values[l] = [
                    i,
                    j,
                    x[j],
                    fixed_coord[1],
                    z[i],
                    array[spacing_i - 1 - i, j],
                ]
                l = l + 1

    if fixed_coord[0] == "Z":
        x = np.linspace(lower_extent[0], upper_extent[0], spacing_j)
        y = np.linspace(lower_extent[1], upper_extent[1], spacing_i)
        l = 0
        for j in range(spacing_j):
            for i in range(spacing_i):
                values[l] = [
                    spacing_i - 1 - i,
                    spacing_j - 1 - j,
                    x[j],
                    y[i],
                    fixed_coord[1],
                    array[spacing_i - 1 - i, j],
                ]
                l = l + 1

    return values
