## Support Functions for intrusion network simulated as the shortest path, and for simulations in general
import numpy as np
import pandas as pd
import random
from LoopStructural.utils import getLogger

logger = getLogger(__name__)

def sort_2_arrays(main_array, array):
    """
    Sort two arrays, considering values of only the main array

    Parameters
    ----------
    main array: numpy array, considered to sort secondary array
    array: numpy aray, array to be sorted

    Returns
    -------
    sorted arrays
    """
    # function to sort 2 arrays, considering values of only the main array

    for i in range(len(main_array)):
        swap = i + np.argmin(main_array[i:])
        (main_array[i], main_array[swap]) = (main_array[swap], main_array[i])
        (array[i], array[swap]) = (array[swap], array[i])

    return main_array, array


def findMinDiff(arr, n):
    """
    Find the min diff by comparing difference of all possible pairs in given array

    Parameters
    ----------
    arr: numpy array with values to compare and fin the minimum difference

    Returns
    -------
    minimum difference between values in arr

    """
    # Initialize difference as infinite
    diff = 10**20

    for i in range(n - 1):
        for j in range(i + 1, n):
            if abs(arr[i] - arr[j]) < diff:
                diff = abs(arr[i] - arr[j])

    return diff


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


def shortest_path(inlet, outlet, time_map):
    """
    Look for the shortest path between inlet and outlet, using a time map.
    In practice, for a given point looks for the neighbour elements which is closer in time.
    The search starts in the inlet, until it reaches the outlet.

    Parameters
    ----------
    inlet: list of indexes (row and column) of inlet point.
    outlet: list of indexes (row and column) of outlet point.
    time_map: numpy array, contains values of time, computed using the fast-marching method.

    Returns
    -------
    inet: (intrusion network) numpy array of same shape of time_map array.
        0 where the intrusion network is found, 1 above it, and -1 below it

    """
    inet = np.ones_like(time_map)  # array to save shortest path with zeros
    temp_inlet = inlet  # temporary inlet
    inet[temp_inlet[0], temp_inlet[1]] = 0
    i = 0

    while True:
        i = i + 1
        time_temp_inlet = time_map[
            temp_inlet[0], temp_inlet[1]
        ]  # obtain time value of temporary outlet
        neighbors = element_neighbour(
            temp_inlet, time_map, inet
        )  # identify neighbours elements of temporary outlet
        direction = index_min(
            neighbors
        )  # obtain the location (index min) of minimun difference
        temp_inlet = new_inlet(temp_inlet, direction)
        a = temp_inlet[0]
        b = temp_inlet[1]

        inet[a, b] = 0

        if temp_inlet[0] == outlet[0] and temp_inlet[1] == outlet[1]:
            break
        else:
            continue

    # Assing -1 to points below intrusion network
    for j in range(len(inet[0])):  # columns
        for h in range(len(inet)):  # rows
            if inet[h, j] == 0:
                break

        inet[(h + 1) :, j] = -1

    return inet


def element_neighbour(index, array, inet):
    """
    Identify the value of neighbours elements for a given element.

    Parameters
    ----------
    index: list of indexes (row and column) of elements.
    array: numpy array, containing values
    inet: numpy array, temporal intrusion network. If one of the elements is already inet=0, the assign -1.

    Returns
    -------
    values: numpy array, 1x5 with the values of the neighbours of a particular element

    """

    rows = len(array) - 1  # max index of rows of time_map array
    cols = len(array[0]) - 1  # max index of columns of time_map arrays
    values = np.zeros(
        9
    )  # array to save values (element above, element to the left, element to the right)
    values[8] = 10
    index_row = index[0]
    index_col = index[1]

    if index_row == 0:
        values[0] = -1
        values[1] = -1
        values[2] = -1

    if index_row == rows:
        values[5] = -1
        values[6] = -1
        values[7] = -1

    if index_col == 0:
        values[0] = -1
        values[3] = -1
        values[5] = -1

    if index_col == cols:
        values[2] = -1
        values[4] = -1
        values[7] = -1

    for k in range(8):
        if values[k] > -1:
            if k == 0:
                values[0] = array[index[0] - 1, index[1] - 1]

            if k == 1:
                values[1] = array[index[0] - 1, index[1]]

            if k == 2:
                values[2] = array[index[0] - 1, index[1] + 1]

            if k == 3:
                values[3] = array[index[0], index[1] - 1]

            if k == 4:
                values[4] = array[index[0], index[1] + 1]

            if k == 5:
                values[5] = array[index[0] + 1, index[1] - 1]

            if k == 6:
                values[6] = array[index[0] + 1, index[1]]

            if k == 7:
                values[7] = array[index[0] + 1, index[1] + 1]

        else:
            continue

    # check if some of the neighbours is already part of the intrusion network
    for h in range(8):
        if values[h] > -1:
            if h == 0:
                if inet[index[0] - 1, index[1] - 1] == 0:
                    values[0] = -1
            if h == 1:
                if inet[index[0] - 1, index[1]] == 0:
                    values[1] = -1
            if h == 2:
                if inet[index[0] - 1, index[1] + 1] == 0:
                    values[2] = -1
            if h == 3:
                if inet[index[0], index[1] - 1] == 0:
                    values[3] = -1
            if h == 4:
                if inet[index[0], index[1] + 1] == 0:
                    values[4] = -1
            if h == 5:
                if inet[index[0] + 1, index[1] - 1] == 0:
                    values[5] = -1
            if h == 6:
                if inet[index[0] + 1, index[1]] == 0:
                    values[6] = -1
            if h == 7:
                if inet[index[0] + 1, index[1] + 1] == 0:
                    values[7] = -1
        else:
            continue

    return values


def index_min(array):
    """
    Given an array of 1x8, dentify the index of the minimum value within the array.

    Parameters
    ----------
    array: numpy array, 1x8

    Returns
    -------
    index_min: integer, index of minimum value in array

    """

    # return the index value of the minimum value in an array of 1x8
    index_array = {}

    for i in range(
        8
    ):  # create a dictionary assining positions from 0 to 7 to the values in the array
        a = i
        if array[i] >= 0:
            index_array.update({i: array[i]})

    minimum_val = min(index_array.values())

    for key, value in index_array.items():
        if value == minimum_val:
            index_min = key

    return index_min


def new_inlet(inlet, direction):
    """
    Determine new inlet indexes, given current inlet and direction of minimum difference in time map.

    Parameters
    ----------
    inlet: list of indexes [row, column]
    direction: integers e[0,7] (0: above-left, 1: above, 2: above right, 3: left, 4: right, 5: below left, 6: below, 7:below right)

    Returns
    -------
    new_inlet: list of indexes [row, column]

    """
    pot_new_inlets = {}

    pot_new_inlets.update({"0": np.array([inlet[0] - 1, inlet[1] - 1])})
    pot_new_inlets.update({"1": np.array([inlet[0] - 1, inlet[1]])})
    pot_new_inlets.update({"2": np.array([inlet[0] - 1, inlet[1] + 1])})
    pot_new_inlets.update({"3": np.array([inlet[0], inlet[1] - 1])})
    pot_new_inlets.update({"4": np.array([inlet[0], inlet[1] + 1])})
    pot_new_inlets.update({"5": np.array([inlet[0] + 1, inlet[1] - 1])})
    pot_new_inlets.update({"6": np.array([inlet[0] + 1, inlet[1]])})
    pot_new_inlets.update({"7": np.array([inlet[0] + 1, inlet[1] + 1])})
    new_outlet = np.zeros(2)

    for key, value in pot_new_inlets.items():
        if key == str(direction):
            new_inlet = value
    return new_inlet


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
