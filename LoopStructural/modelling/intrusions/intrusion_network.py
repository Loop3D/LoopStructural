# import logging
from LoopStructural.utils import getLogger, log_to_file

logger = getLogger(__name__)

try:
    import skfmm as fmm

except ImportError:
    logger.error("skfmm not installed \n" "pip install skfmm")

# from sklearn.cluster import DBSCAN
from LoopStructural.modelling.intrusions.intrusion_support_functions import *

import numpy as np
import pandas as pd
import ckwrap


class IntrusionNetwork:
    def __init__(
        self,
        feature_data=None,
        intrusion_network_contact=None,
        intrusion_network_type=None,
        model=None,
        **kwargs
    ):

        self.feature_data = feature_data
        self.intrusion_network_contact = intrusion_network_contact
        self.intrusion_network_type = intrusion_network_type
        self.model = model
        self.intrusion_network_data = None
        self.other_contact_data = None
        self.anisotropies_series_list = []
        self.anisotropies_series_parameters = {}
        self.anisotropies_fault_list = []
        self.anisotropies_fault_parameters = {}
        self.grid_to_evaluate_ifx = np.zeros([1, 1])
        self.anisotropies_sequence = None
        self.velocity_parameters = None
        self.shortestpath_sections_axis = None
        self.shortestpath_inlet_mean = None
        self.shortestpath_outlet_mean = None
        self.intrusion_network_outcome = None


    def set_data(self):
        """
        Separate data between roof and floor contacts

        Parameters
        ----------
        Returns
        -------
        """
        if self.intrusion_network_contact == "roof":
            other_contact = "floor"
        elif self.intrusion_network_contact == "top":
            other_contact = "base"
        elif self.intrusion_network_contact == "floor":
            other_contact = "roof"
        elif self.intrusion_network_contact == "base":
            other_contact = "top"

        intrusion_network_data = self.feature_data[
            self.feature_data["intrusion_contact_type"]
            == self.intrusion_network_contact
        ]
        other_contact_data = self.feature_data[
            self.feature_data["intrusion_contact_type"] == other_contact
        ]

        self.intrusion_network_data = intrusion_network_data
        self.other_contact_data = other_contact_data

    def set_model(self, model):
        """
        Link a geological model to the feature

        Parameters
        ----------
        model - GeologicalModel

        Returns
        -------

        """
        self.model = model

    def set_contact_anisotropies(self, series_list=None, **kwargs):
        """
        Add to the intrusion network the anisotropies likely exploited by the intrusion (series-type geological features)
        Given a list of series-type features, evaluates contact points on each series and
        compute mean value and standard deviation. Different contacts of the same series are indentify using clustering algortihm.
        Mean and std deviation values will be used to identify each contact thoughout the model using the indicator functions.

        Parameters
        ----------
        series_list = list of series-type features

        Returns
        -------
        assigns to self.anisotropies_series_parameters a list like (for each strtaigraphic contact) =
        [series_name, mean of scalar field vals, standar dev. of scalar field val]

        """
        if self.intrusion_network_type == 'shortest path':
            if "number_contacts" in kwargs:
                n_clusters = kwargs["number_contacts"]
            else:
                n_clusters = [1]*len(series_list)

        if series_list == None:
            self.anisotropies_series_list = None

        else:
            self.anisotropies_series_list = series_list
            series_parameters = {}
            for i in range(len(series_list)):

                data_temp = self.intrusion_network_data[
                    self.intrusion_network_data["intrusion_anisotropy"]
                    == series_list[i].name
                ].copy()
                data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()
                series_i_vals = series_list[i]["feature"].evaluate_value(
                    data_array_temp
                )
                series_array = np.zeros([len(data_array_temp), 4])
                series_array[:, :3] = data_array_temp
                series_array[:, 3] = series_i_vals

                n_contacts = n_clusters[i]

                # use scalar field values to find different contacts
                contact_clustering = ckwrap.ckmeans(series_i_vals, n_contacts)

                for j in range(n_contacts):
                    series_ij_name = series_list[i].name + "_" + str(j)
                    z = np.ma.masked_not_equal(contact_clustering.labels, j)
                    y = np.ma.masked_array(series_i_vals,z.mask)
                    series_ij_vals = np.ma.compressed(y)
                    series_ij_mean = np.mean(series_ij_vals)
                    series_ij_std = np.std(series_ij_vals)
                    
                    series_parameters[series_ij_name] = [
                        series_list[i],
                        series_ij_mean,
                        series_ij_std,
                    ]
                    
            self.anisotropies_series_parameters = series_parameters

    def set_faults_anisotropies(self, fault_list=None):
        """
        Add to the intrusion network the anisotropies likely exploited by the intrusion (fault-type geological features)
        Given a list of fault features, evaluates the contact points on each fault and
        compute mean value and standard deviation.
        These values will be used to identify each fault with the indicator function.

        Parameters
        ----------
        fault_list = list of fault-type features

        Returns
        -------

        """
        if fault_list == None:
            self.anisotropies_fault_list = None
        else:
            self.anisotropies_fault_list = fault_list
            faults_parameters = {}
            for i in range(len(fault_list)):

                data_temp = self.feature_data[self.feature_data["intrusion_anisotropy"]== fault_list[i].name].copy()
                data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()

                fault_i_vals = fault_list[i][0].evaluate_value(data_array_temp)
                fault_i_mean = np.mean(fault_i_vals)
                fault_i_std = np.std(fault_i_vals)

                faults_parameters[fault_list[i].name] = [
                    fault_list[i],
                    fault_i_mean,
                    fault_i_std,
                ]

            self.anisotropies_fault_parameters = faults_parameters

    def create_grid_for_indicator_fxs(self, spacing=None):
        """
        Create the grid points in which to evaluate the indicator functions

        Parameters
        ----------
        spacing = list/array with spacing value for X,Y,Z

        Returns
        -------
        """

        if spacing == None:
            spacing = self.model.nsteps

        grid_points = self.model.regular_grid(spacing, shuffle=False)

        self.grid_to_evaluate_ifx = grid_points

        return grid_points, spacing

    def set_sequence_of_exploited_anisotropies(self, sequence=None):
        """
        Add to the intrusion network the sequence of anisotropies to run the search of the shortest path
        The sequence is used to identify the inlet and outlet for the search of the shortest path.
        If more than two anisotropies, it only considers first and last anisotropies in the sequence.

        Parameters
        ----------

        Returns
        ----------

        """

        self.anisotropies_sequence = sequence

    def set_velocity_parameters(self, velocity_parameters=None):
        """
        Set velocity parameters for the search of the shortest path.
        The velocities should be provided in the same order of the anisotropies sequence.

        Parameters
        ----------

        Returns
        ----------

        """

        self.velocity_parameters = velocity_parameters

    def indicator_function_contacts(self, delta=None):
        """
        Function to compute indicator function for list of contacts anisotropies
        For each point of the grid, this function assignes a 1 if contact i is present, 0 otherwise.

        A contact is defined as an isovalue of an scalar field defining the geological feature  of which the contact is part of.
        Each point is evaluated in the feature scalar field and is identified as part of the contact if its value is around the contact isovalue.

        Parameters
        ----------
        delta : list of numbers, same lenght as number of anisotropies (series). 
            delta multiplies the standard deviation to increase probability of finding the contact on a grid point.

        Returns
        ----------
        Ic: array of [len(grid_points),n_contacts], containing indicator function for list of contacts
        """

        n_series = len(self.anisotropies_series_parameters)  # number of series
        grid_points = self.grid_to_evaluate_ifx

        Ic = np.zeros([len(self.grid_to_evaluate_ifx), n_series])

        delta_list = []

        for i in range(len(delta)):
            for j in range(len(self.anisotropies_series_parameters)):
                delta_list.append(delta[i])


        for i, contact_id in enumerate(self.anisotropies_series_parameters.keys()):
            series_id = self.anisotropies_series_parameters[contact_id][0]
            seriesi_mean = self.anisotropies_series_parameters[contact_id][1]
            seriesi_std = self.anisotropies_series_parameters[contact_id][2]
            seriesi_values = series_id["feature"].evaluate_value(grid_points)

            # apend associated scalar field values to each anisotropy
            self.anisotropies_series_parameters[contact_id].append(seriesi_values)

            # evaluate indicator function in contact (i)
            Ic[np.logical_and(
                (seriesi_mean - seriesi_std * delta_list[i]) <= seriesi_values, 
                seriesi_values <= (seriesi_mean + seriesi_std * delta_list[i])),i] = 1

        return Ic

    def indicator_function_faults(self, delta=None):
        """
        Function to compute indicator function for list of faults anisotropies
        For each point of the grid, this function assignes a 1 if fault i is present, 0 otherwise.

        A fault surface is defined as an isovalue 0 of the scalar field representing the fault surface (coordinate 0 of its structural frame)
        Each point of the grid is evaluated in this scalar field and is identified as part of the fault if its value is around 0.

        Parameters
        ----------
        delta : integer, multiply the standard deviation to increase probability of finding the fault on a point.

        Returns
        ----------
        If: array of [len(grid_points),n_faults], containing indicator function for list of faults
        """

        n_faults = len(self.anisotropies_fault_parameters)  # number of faults
        grid_points = self.grid_to_evaluate_ifx

        If = np.zeros([len(self.grid_to_evaluate_ifx), n_faults])

        for i, fault_id in enumerate(self.anisotropies_fault_parameters.keys()):
            fault_i = self.anisotropies_fault_parameters[fault_id][0]
            faulti_mean = self.anisotropies_fault_parameters[fault_id][1]
            faulti_std = self.anisotropies_fault_parameters[fault_id][2]
            faulti_values = fault_i[0].evaluate_value(grid_points)

            # apend associated scalar field values to each anisotropy
            self.anisotropies_fault_parameters[fault_id].append(faulti_values)
            
            If[np.logical_and(
                (faulti_mean - faulti_std * delta[i]) <= faulti_values,
                faulti_values <= (faulti_mean + faulti_std * delta[i])),i] = 1
            
        return If

    def compute_velocity_field(self, indicator_fx_contacts, indicator_fx_faults):
        """
        Function to compute velocity field to be used in the shortest path search.
        A velocity parameter is assign to each point of the grid, depending on which anisotropy is found at each grid point

        Velocity field
            K(x) = sum(Ic+Kc) + (1-sum(Ic))*sum(If*Kf)
        where Ic and Kc are the indicator function and velocity parameters for contacts,
        and If,Kf are the indicator function and velocity parameters for faults.

        Parameters
        ----------
        spacing = list/array with spacing value for X,Y,Z

        Returns
        -------
        """

        # if no velocity parameters, assign velocities increasing with the order of emplacement
        if self.velocity_parameters == None:
            velocity_parameters = np.arange(len(self.anisotropies_sequence)) + 10
            self.velocity_parameters = velocity_parameters
        else:
            velocity_parameters = self.velocity_parameters.copy()

        # assing velocity parameters k to each anisotropies
        # contacts anisotropies
        IcKc = 0
        IcIc = 0
        anisotropies_sequence_copy = self.anisotropies_sequence.copy()
        if self.anisotropies_series_list != None:
            Kc = np.zeros(len(self.anisotropies_series_parameters))
            for i, key in enumerate(self.anisotropies_series_parameters.keys()):
                series_id = self.anisotropies_series_parameters[key][0]
                for j in range(len(anisotropies_sequence_copy)):
                    if series_id == anisotropies_sequence_copy[j]:
                        Kc[i] = velocity_parameters[j]
                        anisotropies_sequence_copy[j] = None
                        break
            IcKc = np.dot(indicator_fx_contacts, Kc)
            IcIc = np.zeros(len(indicator_fx_contacts))
            for k in range(len(indicator_fx_contacts[0])):
                IcIc = IcIc + indicator_fx_contacts[:, k]

        # faults anisotropies
        IfKf = 0
        if self.anisotropies_fault_list != None:
            Kf = np.zeros(len(self.anisotropies_fault_list))
            for i in range(len(Kf)):
                for j in range(len(velocity_parameters)):
                    if self.anisotropies_fault_list[i] == self.anisotropies_sequence[j]:
                        Kf[i] = velocity_parameters[j]

            IfKf = np.dot(indicator_fx_faults, Kf)

        velocity_field = (IcKc + ((1 - IcIc) * IfKf)) + 0.1

        return velocity_field

    def set_sections_axis(self, axis=None):
        """
        Set section's axis for the search of the shortest path.
        The model volumen in divided in sections parallel to 'axis', and the shortest path is look on each section.

        Parameters
        ----------
        axis: string 'X' or 'Y'

        Returns
        -------
        """

        self.shortestpath_sections_axis = axis

    def build(self, **kwargs):

        # --- check type of intrusion network
        if self.intrusion_network_type == None:
            logger.error(
                "Specify type of intrusion network: 'interpolated' or 'shortest path'"
            )

        elif self.intrusion_network_type == "interpolated":

            inet_points = np.zeros([len(self.intrusion_network_data), 4])
            inet_points_xyz = self.intrusion_network_data.loc[
                :, ["X", "Y", "Z"]
            ].to_numpy()

            inet_points[:, :3] = inet_points_xyz

            self.intrusion_network_outcome = inet_points
            return inet_points

        elif self.intrusion_network_type == "shortest path":

            if "number_contacts" in kwargs:
                n_clusters = kwargs["number_contacts"]
            else:
                n_clusters = 1

            if "delta_c" in kwargs:
                delta_c = kwargs["delta_c"]
            else:
                delta_c = [1] * len(self.anisotropies_series_list)

            if "delta_f" in kwargs:
                delta_f = kwargs["delta_f"]
            else:
                delta_f = [1] * len(self.anisotropies_fault_list)

            grid_points, spacing = self.create_grid_for_indicator_fxs()

            # --- check axis
            if self.shortestpath_sections_axis == "X":
                section_axis = "X"
                other_axis = "Y"
                nxm = spacing[1] * spacing[2]
            elif self.shortestpath_sections_axis == "Y":
                section_axis = "Y"
                other_axis = "X"
                nxm = spacing[0] * spacing[2]
            else:
                logger.error(
                    "Specify axis ('X' or 'Y') to create sections for the search of the shortest path"
                )

            # --- compute indicator functions and velocity field
            Ic = self.indicator_function_contacts(delta=delta_c)

            #--------- check if any anisotropy is indetified by indicator functions:
            if len(np.where(Ic == 1)[0]) == 0:
                logger.error(
                    "No anisotropy identified, increase value of delta_c"
                )

            If = self.indicator_function_faults(delta=delta_f)
            velocity_field = self.compute_velocity_field(Ic, If)

            # --- find first (inlet) and last (outlet) anisotropies in anisotropies sequence, and compute associated scalar fields
            inlet_anisotropy = self.anisotropies_sequence[0]

            if inlet_anisotropy in self.anisotropies_series_list: # if inlet anisotropy type is series
                sf_inlet_anisotropy = inlet_anisotropy["feature"].evaluate_value(grid_points)
   
            else: #otherwise, it is a fault:
                sf_inlet_anisotropy = inlet_anisotropy[0].evaluate_value(grid_points)

            outlet_anisotropy = self.anisotropies_sequence[
                len(self.anisotropies_sequence) - 1
            ]
            if outlet_anisotropy in self.anisotropies_series_list:
                sf_outlet_anisotropy = outlet_anisotropy["feature"].evaluate_value(grid_points)

            else:
                sf_outlet_anisotropy = outlet_anisotropy[0].evaluate_value(grid_points)
            
            # create dataframe containing grid points, velocity values
            nodes = np.linspace(0, len(grid_points), len(grid_points))
            medium = np.zeros([len(grid_points), 7])
            medium[:, 0] = nodes.T
            medium[:, 1:4] = grid_points
            medium[:, 4] = velocity_field
            medium[:, 5] = sf_inlet_anisotropy
            medium[:, 6] = sf_outlet_anisotropy
            medium_df = pd.DataFrame(
                medium,
                columns=["node", "X", "Y", "Z", "k_field", "sf_inlet", "sf_outlet"],
            )

            medium_df.sort_values(
                [other_axis, section_axis, "Z"],
                ascending=[True, False, False],
                inplace=True,
            )

            # Separate data in sections and store them in medium_sections list
            axis_values = medium_df[section_axis].unique()
            medium_sections = []
            for y in range(len(axis_values)):
                section_df = medium_df[medium_df[section_axis] == axis_values[y]].copy()

                section_temp = np.zeros([nxm, 7])
                section_temp[:, 0] = section_df.loc[:, "node"]
                section_temp[:, 1] = section_df.loc[:, "X"]
                section_temp[:, 2] = section_df.loc[:, "Y"]
                section_temp[:, 3] = section_df.loc[:, "Z"]
                section_temp[:, 4] = section_df.loc[:, "k_field"]
                section_temp[:, 5] = section_df.loc[:, "sf_inlet"]
                section_temp[:, 6] = section_df.loc[:, "sf_outlet"]

                medium_sections.append(section_temp)
            
            # Convert each section into arrays containing different set of values
            nodes_arrays = []
            velocity_field_arrays = []
            sf_inlet_arrays = []
            sf_outlet_arrays = []
            inlet_pt_array = []
            outlet_pt_array = []
            initial_point_array = []
            time_maps_array = []
            shortest_path_array = []
            shortest_path_coords = []

            for i in range(len(medium_sections)):
                df_temp = pd.DataFrame(
                    medium_sections[i],
                    columns=["node", "X", "Y", "Z", "k_field", "sf_inlet", "sf_outlet"],
                )
                df_temp.sort_values(
                    [other_axis, "Z"], ascending=[True, False], inplace=True
                )

                # nodes array of each section
                nodes_temp = array_from_coords(df_temp, section_axis, 0)
                nodes_arrays.append(nodes_temp)

                # velocity field arrays
                velocity_field_temp = array_from_coords(df_temp, section_axis, 4)
                velocity_field_arrays.append(velocity_field_temp)

                # array with scalar field associated to first anisotropy (inlet)
                sf_inlet_temp = array_from_coords(df_temp, section_axis, 5)
                sf_inlet_arrays.append(sf_inlet_temp)

                # array with scalar field associated to last anisotropy (outlet)
                sf_outlet_temp = array_from_coords(df_temp, section_axis, 6)
                sf_outlet_arrays.append(sf_outlet_temp)

                inlet_temp, outlet_temp = find_inout_points(
                    velocity_field_temp, self.velocity_parameters
                )
                inlet_pt_array.append(inlet_temp)
                outlet_pt_array.append(outlet_temp)
                initial_pt_array_temp = np.ones(
                    [len(df_temp["Z"].unique()), len(df_temp[other_axis].unique())]
                )
                initial_pt_array_temp[outlet_temp[0], outlet_temp[1]] = 0
                initial_point_array.append(initial_pt_array_temp)

                time_map_temp = fmm.travel_time(
                    initial_pt_array_temp, velocity_field_temp, order=1
                )
                time_maps_array.append(time_map_temp)

                shortest_path_section = shortest_path(
                    inlet_temp, outlet_temp, time_map_temp
                )
                shortest_path_array.append(shortest_path_section)

                shortest_path_coords_section = grid_from_array(
                    shortest_path_section,
                    [section_axis, axis_values[i]],
                    self.model.bounding_box[0],
                    self.model.bounding_box[1],
                )
                shortest_path_coords.append(shortest_path_coords_section)

            # create array of shortest path points
            shortest_path_coords_temp = shortest_path_coords[0]
            for i in range(len(shortest_path_coords) - 1):
                shortest_path_coords_temp = np.vstack(
                    [shortest_path_coords_temp, shortest_path_coords[i + 1]]
                )

            mask = np.ma.masked_not_equal(shortest_path_coords_temp[:,5], 0)
            x = np.ma.compressed(np.ma.masked_array(shortest_path_coords_temp[:,2],mask.mask))
            y = np.ma.compressed(np.ma.masked_array(shortest_path_coords_temp[:,3],mask.mask))
            z = np.ma.compressed(np.ma.masked_array(shortest_path_coords_temp[:,4],mask.mask))
            i = np.ma.compressed(np.ma.masked_array(shortest_path_coords_temp[:,5],mask.mask))
            shortest_path_points = np.array([x,y,z,i]).T
            
            self.intrusion_network_outcome = shortest_path_points
            return shortest_path_points
