# import logging
from LoopStructural.utils import getLogger, log_to_file
logger = getLogger(__name__)
import skfmm as fmm
from LoopStructural.modelling.intrusions.intrusion_support_functions import *

import numpy as np
import pandas as pd
# import random
# from LoopStructural import GeologicalModel


class IntrusionNetwork:
    def __init__(self, feature_data=None, intrusion_network_contact=None, intrusion_network_type = None, model = None,**kwargs):
        
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
        self.grid_to_evaluate_ifx = np.zeros([1,1])
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

        if self.intrusion_network_contact == 'roof':
            other_contact = 'floor'
        elif self.intrusion_network_contact == 'top':
            other_contact = 'base'
        elif self.intrusion_network_contact == 'floor':
            other_contact = 'roof'
        elif self.intrusion_network_contact == 'base':
            other_contact = 'top'

        intrusion_network_data = self.feature_data[self.feature_data['intrusion_contact_type'] == self.intrusion_network_contact]
        other_contact_data = self.feature_data[self.feature_data['intrusion_contact_type'] == other_contact]

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

        
    def set_contact_anisotropies(self, series_list):
        """
        Add to the intrusion network the anisotropies likely exploited by the intrusion (series-type geological features)
        Given a list of series-type features, evaluates the contact points on each series and 
        compute mean value and standard deviation. 
        These values will be used to identify each feature using the indicator functions.

        Parameters
        ----------
        series_list = list of series-type features

        Returns
        -------

        """
        self.anisotropies_series_list = series_list
        series_parameters = {}
        for i in range(len(series_list)):
            sf_series_vals = []
            sf_series_mean = []
            sf_series_std = []
            sf_series_gridpts = []
            
            data_temp = self.intrusion_network_data[self.intrusion_network_data['intrusion_anisotropy'] == series_list[i].name].copy()
            data_array_temp = data_temp.loc[:, ['X','Y','Z']].to_numpy()
            series_i_vals = series_list[i]['feature'].evaluate_value(self.model.scale(data_array_temp, inplace = False))
            
            sf_series_vals.append(series_i_vals)
            sf_series_mean.append(np.mean(series_i_vals))
            sf_series_std.append(np.std(series_i_vals))
            
            series_parameters[series_list[i].name] = [sf_series_vals, sf_series_mean, sf_series_std]
        
        #Pending to add clustering
        
        self.anisotropies_series_parameters = series_parameters
    
    def set_faults_anisotropies(self, fault_list):
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
        self.anisotropies_fault_list = fault_list
        faults_parameters = {}
        for i in range(len(fault_list)):
            sf_faults_vals = []
            sf_faults_mean = []
            sf_faults_std = []
            sf_faults_gridpts = []
            
            data_temp = self.intrusion_network_data[self.intrusion_network_data['intrusion_anisotropy'] == fault_list[i].name].copy()
            data_array_temp = data_temp.loc[:, ['X','Y','Z']].to_numpy()
            fault_i_vals = fault_list[i][0].evaluate_value(self.model.scale(data_array_temp, inplace = False))

            sf_faults_vals.append(fault_i_vals)
            sf_faults_mean.append(np.mean(fault_i_vals))
            sf_faults_std.append(np.std(fault_i_vals))
        
            faults_parameters[fault_list[i].name]  = [sf_faults_vals, sf_faults_mean, sf_faults_std]
        
        self.anisotropies_fault_parameters = faults_parameters
        
        
    def create_grid_for_indicator_fxs(self, spacing = None):
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
            
        x = np.linspace(self.model.origin[0], self.model.maximum[0], spacing[0])
        y = np.linspace(self.model.origin[1], self.model.maximum[1], spacing[1])
        z = np.linspace(self.model.origin[2], self.model.maximum[2], spacing[2])

        xx,yy,zz = np.meshgrid(x,y,z)
        grid_points = np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T
        
        self.grid_to_evaluate_ifx = grid_points
        
        return grid_points, spacing
    
    
    def set_sequence_of_exploited_anisotropies(self, sequence):
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
        
    def set_velocity_parameters(self, velocity_parameters):
        """
        Set velocity parameters for the search of the shortest path. 
        The velocities should be provided in the same order of the anisotropies sequence.

        Parameters
        ----------

        Returns
        ----------

        """
        
        self.velocity_parameters = velocity_parameters

   
    def indicator_function_contacts(self, delta = 1):
        """
        Function to compute indicator function for list of contacts anisotropies
        For each point of the grid, this function assignes a 1 if contact i is present, 0 otherwise.
        
        A contact is defined as an isovalue of an scalar field definin the geological feature  of which the contact is part of.
        Each point is evaluated in the feature scalar field and is identified as part of the contact if its value is around the contact isovalue.

        Parameters
        ----------
        delta : integer, multiply the standard deviation to increase probability of finding the contact on a point.

        Returns
        ----------
        Ic: array of [len(grid_points),n_contacts], containing indicator function for list of contacts
        """

        n_series = len(self.anisotropies_series_parameters)  #number of contacts
        grid_points = self.grid_to_evaluate_ifx
        
        Ic = np.zeros([len(self.grid_to_evaluate_ifx),n_series])
        
        for i in range(n_series):
            series_id = self.anisotropies_series_list[i]
            seriesi_mean = self.anisotropies_series_parameters.get(self.anisotropies_series_list[i].name)[1][0]
            seriesi_std = self.anisotropies_series_parameters.get(self.anisotropies_series_list[i].name)[2][0]
            seriesi_values = series_id['feature'].evaluate_value(self.model.scale(grid_points, inplace = False))
            
            # apend associated scalar field values to each anisotropy
            self.anisotropies_series_parameters.get(self.anisotropies_series_list[i].name).append(seriesi_values)
            
            for j in range(len(seriesi_values)):
                if (seriesi_mean - seriesi_std*delta) <= seriesi_values[j] <= (seriesi_mean + seriesi_std*delta):
                    Ic[j,i] = 1
                else: continue
        
        return Ic
    
    
    def indicator_function_faults(self, delta = 1):
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

        n_faults = len(self.anisotropies_fault_parameters)  #number of faults
        grid_points = self.grid_to_evaluate_ifx
        
        If = np.zeros([len(self.grid_to_evaluate_ifx),n_faults])
        
        for i in range(n_faults):
            fault_id = self.anisotropies_fault_list[i]
            faulti_mean = 0
            faulti_std = 0.005
            faulti_values = fault_id[0].evaluate_value(self.model.scale(grid_points, inplace = False))
            
            # apend associated scalar field values to each anisotropy
            self.anisotropies_fault_parameters.get(self.anisotropies_fault_list[i].name).append(faulti_values)
            
            for j in range(len(faulti_values)):
                if (faulti_mean - faulti_std*delta) <= faulti_values[j] <= (faulti_mean + faulti_std*delta):
                    If[j,i] = 1
                else: continue
        
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
        
        #if no velocity parameters, assign velocities increasing with the order of emplacement
        if self.velocity_parameters == None:
            velocity_parameters = np.zeros(len(self.anisotropies_sequence))
            n=10
            for i in range(len(velocity_parameters)):
                n=n+i
                velocity_parameters[i] = n
            self.velocity_parameters = velocity_parameters
        else:
            velocity_parameters = self.velocity_parameters
            
        # assing velocity parameters k to each anisotropies
        # contacts anisotropies
        IcKc = 0
        IcIc = 0
        if self.anisotropies_series_list != None:
            Kc = np.zeros(len(self.anisotropies_series_list))
            for i in range(len(Kc)):
                for j in range(len(velocity_parameters)):
                    if self.anisotropies_series_list[i] == self.anisotropies_sequence[j]:
                            Kc[i] = velocity_parameters[j]
                            
            IcKc = np.dot(indicator_fx_contacts,Kc)
            IcIc = np.zeros(len(indicator_fx_contacts))
            for k in range(len(indicator_fx_contacts[0])):
                IcIc = IcIc + indicator_fx_contacts[:,k]
            
        
        # faults anisotropies
        IfKf = 0
        if self.anisotropies_fault_list != None:
            Kf = np.zeros(len(self.anisotropies_fault_list))
            for i in range(len(Kf)):
                for j in range(len(velocity_parameters)):
                    if self.anisotropies_fault_list[i] == self.anisotropies_sequence[j]:
                        Kf[i] = velocity_parameters[j]
            
            IfKf = np.dot(indicator_fx_faults,Kf)

        
        velocity_field = (IcKc + ((1-IcIc)*IfKf)) + 0.1     
        
        return velocity_field
    
    def set_sections_axis(self, axis):
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
        
    def build(self, delta_c = 1, delta_f = 1):
        
            
        # --- check type of intrusion network
        if self.intrusion_network_type == None:
            logger.error("Specify type of intrusion network: 'interpolated' or 'shortest path'")
            
        elif self.intrusion_network_type == 'interpolated':
            
            inet_points = np.zeros([len(self.intrusion_network_data),4])
            inet_points[:,:3] = self.intrusion_network_data.loc[:,['X','Y','Z']]
            self.intrusion_network_outcome = inet_points
            return inet_points
        
        elif self.intrusion_network_type == 'shortest path':
            
            grid_points, spacing = self.create_grid_for_indicator_fxs()
            
            # --- check axis
            if self.shortestpath_sections_axis == 'X':
                section_axis = 'X'
                other_axis = 'Y'
                nxm = spacing[1] * spacing[2]
            elif self.shortestpath_sections_axis == 'Y':
                section_axis = 'Y'
                other_axis = 'X'
                nxm = spacing[0] * spacing[2]
            else:
                logger.error("Specify axis ('X' or 'Y') to create sections for the search of the shortest path")
            
            # --- compute indicator functions and velocity field
            Ic = self.indicator_function_contacts(delta = delta_c)
            If = self.indicator_function_faults(delta = delta_f)
            velocity_field = self.compute_velocity_field(Ic,If)
            
            # --- find first (inlet) and last (outlet) anisotropies in anisotropies sequence, and compute associated scalar fields
            inlet_anisotropy = self.anisotropies_sequence[0]
            
            if inlet_anisotropy in self.anisotropies_series_list:
                sf_inlet_anisotropy = inlet_anisotropy['feature'].evaluate_value(self.model.scale(grid_points, inplace = False))
            else:
                sf_inlet_anisotropy = inlet_anisotropy[0].evaluate_value(self.model.scale(grid_points, inplace = False))
            
            outlet_anisotropy = self.anisotropies_sequence[len(self.anisotropies_sequence)-1]
            if outlet_anisotropy in self.anisotropies_series_list:
                sf_outlet_anisotropy = outlet_anisotropy['feature'].evaluate_value(self.model.scale(grid_points, inplace = False))
            else:
                sf_outlet_anisotropy = outlet_anisotropy[0].evaluate_value(self.model.scale(grid_points, inplace = False))
            
            # create dataframe containing grid points, velocity values
            nodes = np.linspace(0,len(grid_points),len(grid_points))
            medium = np.zeros([len(grid_points),7])
            medium[:,0] = nodes.T
            medium[:,1:4] = grid_points
            medium[:,4] = velocity_field
            medium[:,5] = sf_inlet_anisotropy
            medium[:,6] = sf_outlet_anisotropy
            medium_df = pd.DataFrame(medium, columns = ['node','X','Y','Z','k_field','sf_inlet','sf_outlet'])
            
            medium_df.sort_values([other_axis,section_axis,'Z'], ascending = [True, False, False], inplace = True)
            
            # Separate data in sections and store them in medium_sections list
            axis_values = medium_df[section_axis].unique()
            medium_sections = []
            for y in range(len(axis_values)):
                section_df = medium_df[medium_df[section_axis] == axis_values[y]].copy()
                
                section_temp = np.zeros([nxm,7])
                section_temp[:,0] = section_df.loc[: , 'node']
                section_temp[:,1] = section_df.loc[: , 'X']
                section_temp[:,2] = section_df.loc[: , 'Y']
                section_temp[:,3] = section_df.loc[: , 'Z']
                section_temp[:,4] = section_df.loc[: , 'k_field']
                section_temp[:,5] = section_df.loc[: , 'sf_inlet']
                section_temp[:,6] = section_df.loc[: , 'sf_outlet']
                
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
                df_temp = pd.DataFrame(medium_sections[i],columns = ['node','X','Y','Z','k_field','sf_inlet','sf_outlet'])
                df_temp.sort_values([other_axis,'Z'], ascending = [True, False], inplace = True)

                #nodes array of each section
                nodes_temp = array_from_coords(df_temp, section_axis, 0)
                nodes_arrays.append(nodes_temp)
                
                # velocity field arrays 
                velocity_field_temp  = array_from_coords(df_temp, section_axis, 4)
                velocity_field_arrays.append(velocity_field_temp)

                # array with scalar field associated to first anisotropy (inlet)
                sf_inlet_temp = array_from_coords(df_temp, section_axis, 5)
                sf_inlet_arrays.append(sf_inlet_temp)
                
                # array with scalar field associated to last anisotropy (outlet)
                sf_outlet_temp = array_from_coords(df_temp, section_axis, 6)
                sf_outlet_arrays.append(sf_outlet_temp)

                inlet_temp, outlet_temp = find_inout_points(velocity_field_temp, self.velocity_parameters)
                inlet_pt_array.append(inlet_temp)
                outlet_pt_array.append(outlet_temp)
                initial_pt_array_temp = np.ones([len(df_temp['Z'].unique()),len(df_temp[other_axis].unique())])
                initial_pt_array_temp[outlet_temp[0],outlet_temp[1]] = 0
                initial_point_array.append( initial_pt_array_temp)
                
                time_map_temp = fmm.travel_time(initial_pt_array_temp, velocity_field_temp, order = 1)
                time_maps_array.append(time_map_temp)
                
                shortest_path_section = shortest_path(inlet_temp, outlet_temp, time_map_temp)
                shortest_path_array.append(shortest_path_section)
                
                shortest_path_coords_section = grid_from_array(shortest_path_section, [section_axis,axis_values[i]], self.model.origin, self.model.maximum)
                shortest_path_coords.append(shortest_path_coords_section)
            
            #create array of shortest path points
            shortest_path_coords_temp = shortest_path_coords[0]
            for i in range(len(shortest_path_coords)-1):
                shortest_path_coords_temp = np.vstack([shortest_path_coords_temp, shortest_path_coords[i+1]])

            counta = sum(1 for i in range(len(shortest_path_coords_temp)) if shortest_path_coords_temp[i,5] == 0)

            shortest_path_points = np.zeros([counta, 4])
            l = 0
            for k in range(len(shortest_path_coords_temp)):
                if shortest_path_coords_temp[k,5] == 0:
                    shortest_path_points[l,0] = shortest_path_coords_temp[k,2] #X coordinate
                    shortest_path_points[l,1] = shortest_path_coords_temp[k,3] #Y coordinate
                    shortest_path_points[l,2] = shortest_path_coords_temp[k,4] #Z coordinate
                    shortest_path_points[l,3] = shortest_path_coords_temp[k,5] #intrusion network value, must be 0 
                    l=l+1
            self.intrusion_network_outcome = shortest_path_points
            return shortest_path_points

