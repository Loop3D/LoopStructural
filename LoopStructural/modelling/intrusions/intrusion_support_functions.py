## Support Functions for intrusion network simulated as the shortest path, and for simulations in general
import numpy as np
from ...utils import getLogger

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

    if n < 2:
        return diff

    values = np.asarray(arr[:n], dtype=float)
    pairwise_diff = np.abs(values[:, None] - values[None, :])
    np.fill_diagonal(pairwise_diff, np.inf)
    min_diff = pairwise_diff.min()
    if min_diff < diff:
        diff = min_diff

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
        # values are laid out column-major (column j occupies rows
        # n:n+rows of the sorted dataframe, n increasing by rows each column)
        values = df.iloc[:, df_axis].to_numpy()
        array = values.reshape(columns, rows).T
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

    # inlet: leftmost column containing inlet_velocity, take its last (deepest) row match
    inlet_mask = velocity_field_array == inlet_velocity
    col_has_inlet = inlet_mask.any(axis=0)
    if col_has_inlet.any():
        col = int(np.argmax(col_has_inlet))
        rows_matching = np.nonzero(inlet_mask[:, col])[0]
        inlet_point[0] = rows_matching[-1]
        inlet_point[1] = col

    # outlet: rightmost column containing outlet_velocity, take its first row match
    outlet_mask = velocity_field_array == outlet_velocity
    col_has_outlet = outlet_mask.any(axis=0)
    if col_has_outlet.any():
        col = len(col_has_outlet) - 1 - int(np.argmax(col_has_outlet[::-1]))
        rows_matching = np.nonzero(outlet_mask[:, col])[0]
        outlet_point[0] = rows_matching[0]
        outlet_point[1] = col

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

        neighbors = element_neighbour(
            temp_inlet, time_map, inet
        )  # identify neighbours elements of temporary outlet
        direction = index_min(neighbors)  # obtain the location (index min) of minimun difference

        if direction == 10:
            break

        temp_inlet = new_inlet(temp_inlet, direction)
        row = temp_inlet[0]
        col = temp_inlet[1]

        # if row >= n_rows or col >= n_cols:
        #     break
        # else:
        inet[row, col] = 0

        if temp_inlet[0] == outlet[0] and temp_inlet[1] == outlet[1]:
            break
        else:
            continue

    # Assign -1 to points below intrusion network.
    # For each column, find the first row where inet == 0 and set everything
    # below it to -1. Columns with no zero are left untouched (matches the
    # original loop, where h would reach the last row without breaking and
    # inet[(h + 1):, j] = -1 is then a no-op empty slice).
    mask_zero = inet == 0
    has_zero = mask_zero.any(axis=0)
    first_zero_row = np.argmax(mask_zero, axis=0)
    row_idx = np.arange(inet.shape[0])[:, None]
    below_mask = (row_idx > first_zero_row[None, :]) & has_zero[None, :]
    inet[below_mask] = -1

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

    # fixed offsets of the 8 neighbours, in the same order as the original
    # k/h indices (0: above-left, 1: above, 2: above-right, 3: left, 4: right,
    # 5: below-left, 6: below, 7: below-right)
    offsets = np.array(
        [
            [-1, -1],
            [-1, 0],
            [-1, 1],
            [0, -1],
            [0, 1],
            [1, -1],
            [1, 0],
            [1, 1],
        ]
    )
    neighbour_idx = np.asarray(index) + offsets
    valid = (
        (neighbour_idx[:, 0] >= 0)
        & (neighbour_idx[:, 0] <= rows)
        & (neighbour_idx[:, 1] >= 0)
        & (neighbour_idx[:, 1] <= cols)
    )

    values = np.full(8, -1.0)
    if valid.any():
        valid_rows = neighbour_idx[valid, 0]
        valid_cols = neighbour_idx[valid, 1]
        values[valid] = array[valid_rows, valid_cols]

        # check if some of the neighbours is already part of the intrusion network
        already_in_network = inet[valid_rows, valid_cols] == 0
        valid_positions = np.nonzero(valid)[0]
        values[valid_positions[already_in_network]] = -2

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
    array = np.asarray(array)
    mask = array >= 0

    if mask.any():
        masked = np.where(mask, array, np.inf)
        minimum_val = masked.min()
        # original loop keeps overwriting index_min for every matching key
        # in increasing order, so ties resolve to the LAST (highest) index
        matches = np.nonzero(masked == minimum_val)[0]
        index_min = int(matches[-1])
    else:
        index_min = 10

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
    new_inlet = np.zeros(2)

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

    array = np.asarray(array)
    spacing_i = len(array)  # number of rows
    spacing_j = len(array[0])  # number of columns
    values = np.zeros([spacing_i * spacing_j, 6])

    # original loops iterate outer j, inner i, with l incrementing each
    # inner step, so i is the fast-varying axis and j the slow-varying axis
    i_flat = np.tile(np.arange(spacing_i), spacing_j)
    j_flat = np.repeat(np.arange(spacing_j), spacing_i)
    array_vals = array[spacing_i - 1 - i_flat, j_flat]

    if fixed_coord[0] == "X":
        y = np.linspace(lower_extent[1], upper_extent[1], spacing_j)
        z = np.linspace(lower_extent[2], upper_extent[2], spacing_i)
        values[:, 0] = i_flat
        values[:, 1] = j_flat
        values[:, 2] = fixed_coord[1]
        values[:, 3] = y[j_flat]
        values[:, 4] = z[i_flat]
        values[:, 5] = array_vals

    if fixed_coord[0] == "Y":
        x = np.linspace(lower_extent[0], upper_extent[0], spacing_j)
        z = np.linspace(lower_extent[2], upper_extent[2], spacing_i)
        values[:, 0] = i_flat
        values[:, 1] = j_flat
        values[:, 2] = x[j_flat]
        values[:, 3] = fixed_coord[1]
        values[:, 4] = z[i_flat]
        values[:, 5] = array_vals

    if fixed_coord[0] == "Z":
        x = np.linspace(lower_extent[0], upper_extent[0], spacing_j)
        y = np.linspace(lower_extent[1], upper_extent[1], spacing_i)
        values[:, 0] = spacing_i - 1 - i_flat
        values[:, 1] = spacing_j - 1 - j_flat
        values[:, 2] = x[j_flat]
        values[:, 3] = y[i_flat]
        values[:, 4] = fixed_coord[1]
        values[:, 5] = array_vals

    return values
