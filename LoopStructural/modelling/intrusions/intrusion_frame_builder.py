from ...modelling.features.builders import StructuralFrameBuilder
from ...modelling.features import BaseFeature
from ...modelling.intrusions.intrusion_support_functions import (
    grid_from_array,
    shortest_path,
    array_from_coords,
    find_inout_points,
)
from ...utils import getLogger


logger = getLogger(__name__)

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

try:
    import skfmm as fmm

except ImportError:
    logger.warning(
        "Cannot use IntrusionFrameBuilder: skfmm not installed \n"
        "pip install scikit-fmm"
    )


class IntrusionFrameBuilder(StructuralFrameBuilder):
    def __init__(self, interpolator=None, interpolators=None, model=None, **kwargs):
        """IntrusionBuilder set up the intrusion frame to build an intrusion
            The intrusion frame is curvilinear coordinate system of the intrusion that controls the simulation of the intrusion extent.
            The object is constrained with intrusion network points (computed) and flow and inflation measurements (provided by the user).
            The intrusion network is a representation of the approximated location of roof (or floor) contact of the intrusion.
            The intrusion network be constrained using the anisotropies of the host rock if the roof (or floor) contact is not well constrained.

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        model : GeologicalModel
            reference to the model containing the fault
        """

        StructuralFrameBuilder.__init__(self, interpolator, interpolators, **kwargs)

        self.origin = np.array([np.nan, np.nan, np.nan])
        self.maximum = np.array([np.nan, np.nan, np.nan])
        self.model = model
        self.minimum_origin = self.model.bounding_box[0, :]
        self.maximum_maximum = self.model.bounding_box[1, :]

        # -- intrusion network input data
        self.intrusion_network_contact = None
        self.intrusion_network_type = None
        self.intrusion_network_data = None
        self.other_contact_data = None
        self.post_intrusion_faults = None
        self.grid_to_evaluate_ifx = np.zeros([1, 1])

        self.anisotropies_series_list = []
        self.anisotropies_series_parameters = {}
        self.anisotropies_fault_list = []
        self.anisotropies_fault_parameters = {}

        self.anisotropies_sequence = None
        self.velocity_parameters = None
        self.shortestpath_sections_axis = None
        self.number_of_contacts = None
        self.delta_contacts = None
        self.delta_faults = None
        # self.steps_stratigraphic_units = None
        self.intrusion_steps = None
        self.intrusion_network_points = None
        self.intrusion_network_gradients = None

        self.velocity_field_arrays = None
        self.IFf = None
        self.IFc = None

    def update_geometry(self, points):
        self.origin = np.nanmin(np.array([np.min(points, axis=0), self.origin]), axis=0)
        self.maximum = np.nanmax(
            np.array([np.max(points, axis=0), self.maximum]), axis=0
        )
        self.origin[self.origin < self.minimum_origin] = self.minimum_origin[
            self.origin < self.minimum_origin
        ]
        self.maximum[self.maximum > self.maximum_maximum] = self.maximum_maximum[
            self.maximum > self.maximum_maximum
        ]

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

    def set_intrusion_network_data(self, intrusion_data):
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

        intrusion_network_data = intrusion_data[
            intrusion_data["intrusion_contact_type"] == self.intrusion_network_contact
        ]
        other_contact_data = intrusion_data[
            intrusion_data["intrusion_contact_type"] == other_contact
        ]

        self.intrusion_network_data = intrusion_network_data
        self.other_contact_data = other_contact_data

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

    def add_contact_anisotropies(self, series_list=[], **kwargs):
        """
        Add to the intrusion network the anisotropies likely
        exploited by the intrusion (series-type geological features)
        Given a list of series-type features, evaluates contact points
        on each series and compute mean value and standard deviation.
        Different contacts of the same series are indentify using
        clustering algortihm. Mean and std deviation values will
        be used to identify each contact thoughout the model using
        the indicator functions.

        Parameters
        ----------
        series_list = list of series-type features

        Returns
        -------

        Note
        -----
        assigns to self.anisotropies_series_parameters a list like (for each strtaigraphic contact) =
        [series_name, mean of scalar field vals, standar dev. of scalar field val]

        """
        if self.intrusion_network_type == "shortest path":
            n_clusters = self.number_of_contacts

        self.anisotropies_series_list = series_list
        series_parameters = {}
        for i, series in enumerate(series_list):

            series.faults_enabled = True  # data points evaluated in faulted hostrock
            data_temp = self.intrusion_network_data[
                self.intrusion_network_data["intrusion_anisotropy"] == series.name
            ].copy()
            data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()
            series_i_vals = series.evaluate_value(data_array_temp)
            series_array = np.zeros((len(data_array_temp), 4))
            series_array[:, :3] = data_array_temp
            series_array[:, 3] = series_i_vals

            n_contacts = n_clusters[i]

            # use scalar field values to find different contacts
            series_i_vals_mod = series_i_vals.reshape(len(series_i_vals), 1)
            # TODO create global loopstructural random state variable
            contact_clustering = KMeans(n_clusters=n_contacts, random_state=0).fit(
                series_i_vals_mod
            )

            for j in range(n_contacts):
                z = np.ma.masked_not_equal(contact_clustering.labels_, j)
                y = np.ma.masked_array(series_i_vals, z.mask)
                series_ij_vals = np.ma.compressed(y)
                series_ij_mean = np.mean(series_ij_vals)
                series_ij_std = np.std(series_ij_vals)
                series_ij_name = f"{series.name}_{str(series_ij_mean)}"

                series_parameters[series_ij_name] = [
                    series,
                    series_ij_mean,
                    series_ij_std,
                ]

            self.anisotropies_series_parameters = series_parameters

    def add_faults_anisotropies(self, fault_list=[]):
        """
        Add to the intrusion network the anisotropies likely
        exploited by the intrusion (fault-type geological features)
        Given a list of fault features, evaluates the contact
        points on each fault and compute mean value and standard
        deviation. These values will be used to identify each
        fault with the indicator function.

        Parameters
        ----------
        fault_list = list of fault-type features

        Returns
        -------

        """
        self.anisotropis_fault_list = fault_list

        faults_parameters = {}
        for fault in fault_list:

            data_temp = self.intrusion_network_data[
                self.intrusion_network_data["intrusion_anisotropy"] == fault.name
            ].copy()
            data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()

            if data_temp.empty:
                fault_i_mean = 0
                fault_i_std = 0.1

            else:
                # evaluate the value of the fault, evaluate_value called on a fault
                # will return the scalar field value of the main structural feature
                fault_i_vals = fault.evaluate_value(data_array_temp)
                fault_i_mean = np.mean(fault_i_vals)
                fault_i_std = np.std(fault_i_vals)

            faults_parameters[fault.name] = [
                fault,
                fault_i_mean,
                fault_i_std,
            ]

        self.anisotropies_fault_parameters = faults_parameters

    def add_steps_stratigraphic_units(self):
        """
        adds stratigraphic contacts from before and after a step (controlled by a fault)
        for each fault, at least one strigraphic unit must be provided.


        Parameters
        ----------
        fault_list = list of fault-type features

        Returns
        -------

        """

        if self.model.stratigraphic_column is None:
            logger.error(
                "No stratigraphic column. \
                Set stratigraphic column using \
                model.set_stratigraphic_column(name_of_stratigraphic_column)"
            )

        # set data
        intrusion_network_data = self.intrusion_network_data.copy()
        intrusion_network_data_xyz = (
            intrusion_network_data.loc[:, ["X", "Y", "Z"]].copy().to_numpy()
        )
        std_backup = 25

        for step_i, step in self.intrusion_steps.items():
            step_structure = step.get("structure")
            unit_from_name = step.get("unit_from")
            series_from_name = step.get("series_from")
            unit_from_id = self.model.stratigraphic_column[series_from_name.name][
                unit_from_name
            ].get("id")

            unit_to_name = step.get("unit_to")
            series_to_name = step.get("series_to")
            unit_to_id = self.model.stratigraphic_column[series_to_name.name][
                unit_to_name
            ].get("id")

            series_from_name.faults_enabled = (
                True  # data points evaluated in faulted hostrock
            )
            series_to_name.faults_enabled = (
                True  # data points evaluated in faulted hostrock
            )

            intrusion_network_data.loc[:, "model_values"] = self.model.evaluate_model(
                self.model.rescale(intrusion_network_data_xyz, inplace=False)
            )

            # check if step is within the same unit:

            if unit_from_name == unit_to_name:
                data_points_xyz = (
                    intrusion_network_data[
                        intrusion_network_data["model_values"] == unit_from_id
                    ]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                series_values = series_from_name.evaluate_value(data_points_xyz)
                series_values_mod = series_values.reshape(len(series_values), 1)
                contact_clustering = KMeans(n_clusters=2, random_state=0).fit(
                    series_values_mod
                )

                # contact 0
                z = np.ma.masked_not_equal(contact_clustering.labels_, 0)
                y = np.ma.masked_array(series_values, z.mask)
                contact_0_vals = np.ma.compressed(y)
                contact_0_mean = np.mean(contact_0_vals)
                contact_0_std = np.std(contact_0_vals)

                # contact 1
                z = np.ma.masked_not_equal(contact_clustering.labels_, 1)
                y = np.ma.masked_array(series_values, z.mask)
                contact_1_vals = np.ma.compressed(y)
                contact_1_mean = np.mean(contact_1_vals)
                contact_1_std = np.std(contact_1_vals)

                if contact_0_mean <= contact_1_mean:
                    step["unit_from_mean"] = contact_0_mean
                    step["unit_from_std"] = contact_0_std
                    step["unit_to_mean"] = contact_1_mean
                    step["unit_to_std"] = contact_1_std

                else:
                    step["unit_from_mean"] = contact_1_mean
                    step["unit_from_std"] = contact_1_std
                    step["unit_to_mean"] = contact_0_mean
                    step["unit_to_std"] = contact_0_std

            else:  # step between different strat units
                data_points_from_xyz = (
                    intrusion_network_data[
                        intrusion_network_data["model_values"] == unit_from_id
                    ]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                step_structure_points_vals = step_structure[0].evaluate_value(
                    data_points_from_xyz
                )
                # print(step_structure_points_vals)
                if len(data_points_from_xyz) == 0:  # no data points in strat unit
                    unit_from_min = self.model.stratigraphic_column[
                        series_from_name.name
                    ][unit_from_name].get("min")
                    unit_from_max = self.model.stratigraphic_column[
                        series_from_name.name
                    ][unit_from_name].get("max")
                    step["unit_from_mean"] = (
                        unit_from_min + (unit_from_max - unit_from_min) / 2
                    )
                    step["unit_from_std"] = std_backup
                else:
                    series_values = series_from_name.evaluate_value(
                        data_points_from_xyz
                    )
                    # print(len(step_structure_points_vals), len(series_values))
                    mask = step_structure_points_vals < 0
                    if len(mask) > 0:
                        series_values_mod = np.ma.compressed(
                            np.ma.masked_array(
                                series_values, step_structure_points_vals < 0
                            )
                        )
                    else:
                        series_values_mod = series_values
                    step["unit_from_mean"] = np.nanmean(series_values_mod)
                    step["unit_from_std"] = np.nanstd(series_values_mod)

                if step["unit_from_std"] == 0:
                    step["unit_from_std"] = std_backup

                data_points_to_xyz = (
                    intrusion_network_data[
                        intrusion_network_data["model_values"] == unit_to_id
                    ]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                step_structure_points_vals = step_structure[0].evaluate_value(
                    data_points_to_xyz
                )
                if len(data_points_to_xyz) == 0:
                    unit_to_min = self.model.stratigraphic_column[series_to_name.name][
                        unit_to_name
                    ].get("min")
                    unit_to_max = self.model.stratigraphic_column[series_to_name.name][
                        unit_to_name
                    ].get("max")
                    step["unit_to_mean"] = unit_to_min + (unit_to_max - unit_to_min) / 2
                    step["unit_to_std"] = std_backup
                else:
                    series_values = series_to_name.evaluate_value(data_points_to_xyz)
                    # print(len(step_structure_points_vals), len(series_values))
                    mask = step_structure_points_vals > 0
                    if len(mask) > 0:
                        series_values_mod = np.ma.compressed(
                            np.ma.masked_array(
                                series_values, step_structure_points_vals > 0
                            )
                        )
                    else:
                        series_values_mod = series_values
                    step["unit_to_mean"] = np.nanmean(series_values_mod)
                    step["unit_to_std"] = np.nanstd(series_values_mod)

                if step["unit_to_std"] == 0:
                    step["unit_to_std"] = std_backup

    def set_intrusion_network_parameters(
        self, intrusion_data, intrusion_network_input, **kwargs
    ):

        """
        Set variables to create intrusion network.

        Parameters
        ----------
        intrusion_data = DataFrame, intrusion contact data
        intrusion_network_input = Dictionary,
            contact : string, contact of the intrusion to be used to create the network (roof or floor)
            type : string, type of algorithm to create the intrusion network (interpolated or shortest path).
                    Shortest path is recommended when intrusion contact is not well constrained
            contacts_anisotropies : list of series-type features involved in intrusion emplacement
            structures_anisotropies : list of fault-type features involved in intrusion emplacement
            sequence_anisotropies : list of anisotropies to look for the shortest path. It could be only starting and end point.

        Returns
        -------

        """

        self.intrusion_network_contact = intrusion_network_input.get("contact", "floor")

        self.intrusion_network_type = intrusion_network_input.get(
            "type", "interpolated"
        )

        self.set_intrusion_network_data(intrusion_data)

        if self.intrusion_network_type == "interpolated":

            intrusion_steps = intrusion_network_input.get("intrusion_steps", None)
            if intrusion_steps is not None:

                self.intrusion_steps = intrusion_steps

                self.delta_contacts = intrusion_network_input.get(
                    "delta_c", [1] * len(intrusion_steps)
                )

                self.delta_faults = intrusion_network_input.get(
                    "delta_f", [1] * len(intrusion_steps)
                )

                self.add_steps_stratigraphic_units()  # update steps parameters

                # set fault anisotropies and compute parameters
                fault_anisotropies = intrusion_network_input.get(
                    "structures_anisotropies", None
                )

                if fault_anisotropies is None:
                    fault_anisotropies = []
                    for step in self.intrusion_steps.keys():
                        fault_anisotropies.append(
                            self.intrusion_steps[step].get("structure")
                        )

                    self.add_faults_anisotropies(fault_anisotropies)

            # add contact anisotropies list
            contact_anisotropies = intrusion_network_input.get(
                "contact_anisotropies", None
            )
            self.anisotropies_series_list = contact_anisotropies

        if self.intrusion_network_type == "shortest path":

            # set the sequence of anisotropies to follow by the shortest path algorithm
            self.anisotropies_sequence = intrusion_network_input.get(
                "shortest_path_sequence", []
            )

            # set velocity parameters for the search of the shortest path
            self.velocity_parameters = intrusion_network_input.get(
                "velocity_parameters", []
            )

            # set axis for shortest path algorithm
            self.shortestpath_sections_axis = intrusion_network_input.get(
                "shortest_path_axis", []
            )

            # add contact anisotropies and compute parameters for shortest path algorithm
            contact_anisotropies = intrusion_network_input.get(
                "contact_anisotropies", []
            )
            # set number of contacts
            self.number_of_contacts = intrusion_network_input.get(
                "number_contacts", [1] * len(contact_anisotropies)
            )
            self.add_contact_anisotropies(contact_anisotropies)

            # set fault anisotropies and compute parameters for shortest path algorithm
            fault_anisotropies = intrusion_network_input.get(
                "structures_anisotropies", []
            )
            self.add_faults_anisotropies(fault_anisotropies)

            # set delta_c for indicator function
            self.delta_contacts = intrusion_network_input.get(
                "delta_c", [1] * np.sum(np.array(self.number_of_contacts))
            )

            # set delta_f for indicator function
            self.delta_faults = intrusion_network_input.get(
                "delta_f", [1] * len(self.anisotropies_fault_list)
            )

    def indicator_function_contacts(self, delta=None):
        """
        Function to compute indicator function for list of contacts anisotropies
        For each point of the grid, this function assignes a 1 if contact i is present, 0 otherwise.

        A contact is defined as an isovalue of an scalar field defining the
        geological feature  of which the contact is part of.
        Each point is evaluated in the feature scalar field and is
        identified as part of the contact if its value is around the contact isovalue.

        Parameters
        ----------
        delta : list of numbers, same lenght as number of anisotropies (series).
            delta multiplies the standard deviation to increase
            probability of finding the contact on a grid point.

        Returns
        ----------
        Ic: array of [len(grid_points),n_contacts], containing indicator function for list of contacts
        """

        n_series = len(self.anisotropies_series_parameters)  # number of series
        grid_points = self.grid_to_evaluate_ifx

        Ic = np.zeros([len(self.grid_to_evaluate_ifx), n_series])

        delta_list = delta

        for i, contact_id in enumerate(sorted(self.anisotropies_series_parameters)):
            series_id = self.anisotropies_series_parameters[contact_id][0]
            seriesi_mean = self.anisotropies_series_parameters[contact_id][1]
            seriesi_std = self.anisotropies_series_parameters[contact_id][2]

            series_id.faults_enabled = True
            seriesi_values = series_id.evaluate_value(grid_points)

            # apend associated scalar field values to each anisotropy
            self.anisotropies_series_parameters[contact_id].append(seriesi_values)

            # evaluate indicator function in contact (i)
            Ic[
                np.logical_and(
                    (seriesi_mean - seriesi_std * delta_list[i]) <= seriesi_values,
                    seriesi_values <= (seriesi_mean + seriesi_std * delta_list[i]),
                ),
                i,
            ] = 1

        self.IFc = Ic
        return Ic

    def indicator_function_faults(self, delta=None):
        """
        Function to compute indicator function for list of faults anisotropies
        For each point of the grid, this function assignes a 1 if fault i is present, 0 otherwise.

        A fault surface is defined as an isovalue 0 of the scalar field representing
         the fault surface (coordinate 0 of its structural frame)
        Each point of the grid is evaluated in this scalar field
        and is identified as part of the fault if its value is around 0.

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

            If[
                np.logical_and(
                    (faulti_mean - faulti_std * delta[i]) <= faulti_values,
                    faulti_values <= (faulti_mean + faulti_std * delta[i]),
                ),
                i,
            ] = 1

        self.IFf = If
        return If

    def compute_velocity_field(self, indicator_fx_contacts, indicator_fx_faults):
        """
        Function to compute velocity field to be used in the shortest path search.
        A velocity parameter is assign to each point of the grid, depending on
        which anisotropy is found at each grid point

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
        if self.velocity_parameters == None or len(self.velocity_parameters) == 0:
            velocity_parameters = (
                np.arange(len(self.anisotropies_sequence)) + 10
            )  # TODO why +10?
            self.velocity_parameters = velocity_parameters
        else:
            velocity_parameters = self.velocity_parameters.copy()

        # assing velocity parameters k to each anisotropies
        # contacts anisotropies
        IcKc = 0
        IcIc = 0
        anisotropies_sequence_copy = self.anisotropies_sequence.copy()
        if self.anisotropies_series_list:
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

    def create_intrusion_network(self, **kwargs) -> np.ndarray:

        """
        Created a numpy array containing (x,y,z) coordinates of intrusion network points

        Parameters
        ----------

        Returns
        -------
        intrusion_network_points : np.ndarray
        """

        # --- check type of intrusion network
        if not self.intrusion_network_type:
            logger.error(
                "Specify type of intrusion network: 'interpolated' or 'shortest path'"
            )

        elif self.intrusion_network_type == "interpolated":

            # 1. data points of intrusion network contacts
            # intrusion_network_points = np.zeros([len(self.intrusion_network_data), 4])
            inet_points_xyz = self.intrusion_network_data.loc[
                :, ["X", "Y", "Z"]
            ].to_numpy()

            intrusion_network_points = inet_points_xyz

            grid_points, spacing = self.create_grid_for_indicator_fxs()

            # self.intrusion_network_points = intrusion_network_points

            # 2. More constraints if steps present:
            if self.intrusion_steps is not None:

                If = self.indicator_function_faults(delta=self.delta_faults)

                if len(np.where(If == 1)[0]) == 0:
                    logger.error("No faults identified, increase value of delta_f")

                If_sum = np.sum(If, axis=1)  # .reshape(len(If),1)

                # evaluate grid points in series

                for i, step_i in enumerate(self.intrusion_steps.keys()):
                    delta_contact = self.delta_contacts[i]
                    if type(delta_contact) is list:
                        delta_contact0 = delta_contact[0]
                        delta_contact1 = delta_contact[1]
                    else:
                        delta_contact0 = delta_contact
                        delta_contact1 = delta_contact

                    step_fault = self.intrusion_steps[step_i].get("structure")
                    series_from_name = self.intrusion_steps[step_i].get("series_from")
                    series_to_name = self.intrusion_steps[step_i].get("series_to")

                    # TODO LG to ask FA why this is done? aren't faults
                    # enabled by default? If they ar disabled, should
                    # they get turned back off after this?
                    series_from_name.faults_enabled = True
                    series_to_name.faults_enabled = True

                    fault_gridpoints_vals = step_fault[0].evaluate_value(grid_points)

                    series_from_gridpoints_vals = series_from_name.evaluate_value(
                        grid_points
                    )
                    # if
                    if series_from_name == series_to_name:
                        series_to_gridpoints_vals = series_from_gridpoints_vals

                    else:
                        series_to_gridpoints_vals = series_to_name.evaluate_value(
                            grid_points
                        )

                    contacts0_val_min = self.intrusion_steps[step_i].get(
                        "unit_from_mean"
                    ) - (
                        self.intrusion_steps[step_i].get("unit_from_std")
                        * delta_contact0
                    )
                    contacts0_val_max = self.intrusion_steps[step_i].get(
                        "unit_from_mean"
                    ) + (
                        self.intrusion_steps[step_i].get("unit_from_std")
                        * delta_contact0
                    )
                    contacts1_val_min = self.intrusion_steps[step_i].get(
                        "unit_to_mean"
                    ) - (
                        self.intrusion_steps[step_i].get("unit_to_std") * delta_contact1
                    )
                    contacts1_val_max = self.intrusion_steps[step_i].get(
                        "unit_to_mean"
                    ) + (
                        self.intrusion_steps[step_i].get("unit_to_std") * delta_contact1
                    )

                    self.intrusion_steps[step_i]["constraints_hw"] = grid_points[
                        np.logical_and(
                            If_sum == 0,
                            (fault_gridpoints_vals >= 0)
                            & (series_from_gridpoints_vals >= contacts0_val_min)
                            & (series_from_gridpoints_vals <= contacts0_val_max),
                        )
                    ]

                    self.intrusion_steps[step_i]["constraints_fw"] = grid_points[
                        np.logical_and(
                            If_sum == 0,
                            (fault_gridpoints_vals < 0)
                            & (series_from_gridpoints_vals >= contacts0_val_min)
                            & (series_from_gridpoints_vals <= contacts0_val_max),
                        )
                    ]

                    fault_i_points = grid_points[
                        np.logical_and(
                            If_sum == 0,
                            np.logical_or(
                                (fault_gridpoints_vals >= 0)
                                & (series_from_gridpoints_vals >= contacts0_val_min)
                                & (series_from_gridpoints_vals <= contacts0_val_max),
                                (fault_gridpoints_vals < 0)
                                & (series_to_gridpoints_vals >= contacts1_val_min)
                                & (series_to_gridpoints_vals <= contacts1_val_max),
                            ),
                        )
                    ]

                    # #TODO
                    # gradients_in_fault =

                    intrusion_network_points = np.vstack(
                        [intrusion_network_points, fault_i_points]
                    )

                    splits_from_sill_name = self.intrusion_steps[step_i].get(
                        "splits_from", None
                    )

                    # check if sill comes from anothe sill
                    # (add all the constraint from original sill, this have to be changed and adapted so it adds constraints for specific faults)

                    if splits_from_sill_name is not None:
                        splits_from_sill_steps = self.model.__getitem__(
                            splits_from_sill_name
                        ).intrusion_frame.builder.intrusion_steps
                        for j, step_j in enumerate(splits_from_sill_steps.keys()):
                            step_j_hg_constraints = splits_from_sill_steps[step_j].get(
                                "constraints_hw"
                            )
                            intrusion_network_points = np.vstack(
                                [intrusion_network_points, step_j_hg_constraints]
                            )

            self.intrusion_network_points = intrusion_network_points

            # 3. Gradient contraints for intrusion frame coord0
            # Evaluate points in gradient of stratigraphy, and exclude points around faults
            if self.anisotropies_series_list:
                series_id = self.anisotropies_series_list[0]
                if not isinstance(series_id, BaseFeature):
                    raise TypeError(
                        f"Cannot set anisotropy with {type(series_id)}, must be a geological feature"
                    )
                series_id.faults_enabled = True
                strat_gradient_grid_points = series_id.evaluate_gradient(grid_points)

                # If intrusion network is built usinf roof/top contact, then change vector direction
                if (
                    self.intrusion_network_contact == "roof"
                    or self.intrusion_network_contact == "top"
                ):
                    grid_points_inflation = strat_gradient_grid_points * (-1)
                else:
                    grid_points_inflation = strat_gradient_grid_points

                grid_points_and_inflation_all = np.hstack(
                    [grid_points, grid_points_inflation]
                )

                if self.intrusion_steps is None:
                    If_sum = np.zeros(
                        len(grid_points_and_inflation_all)
                    ).T  # e.g. no faults afecting the intrusion

                grid_points_and_inflation = grid_points_and_inflation_all[If_sum == 0]
                self.intrusion_network_gradients = grid_points_and_inflation

            return intrusion_network_points

        elif self.intrusion_network_type == "shortest path":

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
            Ic = self.indicator_function_contacts(delta=self.delta_contacts)

            # --------- check if any anisotropy is indetified by indicator functions:
            if len(np.where(Ic == 1)[0]) == 0:
                logger.error("No anisotropy identified, increase value of delta_c")

            If = self.indicator_function_faults(delta=self.delta_faults)
            # --------- check if any fault is indetified by indicator functions:
            if len(np.where(If == 1)[0]) == 0:
                logger.error("No faults identified, increase value of delta_f")

            velocity_field = self.compute_velocity_field(Ic, If)

            # --- find first (inlet) and last (outlet) anisotropies in anisotropies sequence, and compute associated scalar fields
            inlet_anisotropy = self.anisotropies_sequence[0]

            if (
                inlet_anisotropy in self.anisotropies_series_list
            ):  # if inlet anisotropy type is series
                sf_inlet_anisotropy = inlet_anisotropy.evaluate_value(grid_points)

            else:  # otherwise, it is a fault:
                sf_inlet_anisotropy = inlet_anisotropy[0].evaluate_value(grid_points)

            outlet_anisotropy = self.anisotropies_sequence[
                len(self.anisotropies_sequence) - 1
            ]
            if outlet_anisotropy in self.anisotropies_series_list:
                sf_outlet_anisotropy = outlet_anisotropy.evaluate_value(grid_points)

            else:
                sf_outlet_anisotropy = outlet_anisotropy[0].evaluate_value(grid_points)

            # create dataframe containing grid points, velocity values
            nodes = np.linspace(0, len(grid_points), len(grid_points))
            medium = np.zeros((len(grid_points), 7))
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

            mask = np.ma.masked_not_equal(shortest_path_coords_temp[:, 5], 0)
            x = np.ma.compressed(
                np.ma.masked_array(shortest_path_coords_temp[:, 2], mask.mask)
            )
            y = np.ma.compressed(
                np.ma.masked_array(shortest_path_coords_temp[:, 3], mask.mask)
            )
            z = np.ma.compressed(
                np.ma.masked_array(shortest_path_coords_temp[:, 4], mask.mask)
            )
            i = np.ma.compressed(
                np.ma.masked_array(shortest_path_coords_temp[:, 5], mask.mask)
            )
            shortest_path_points = np.array([x, y, z, i]).T

            self.intrusion_network_points = shortest_path_points

            self.velocity_field_arrays = velocity_field_arrays
            return shortest_path_points

    def get_indicator_function_points(self, ifx_type="contacts"):
        """returns the points indicated as p[art of contact or fault anisotropies.

        Parameters
        ----------
        ifx_type : string,
                    'contacts' or 'faults'

        """

        grid_points = self.grid_to_evaluate_ifx

        if ifx_type == "contacts":
            IF = self.IFc
        else:
            IF = self.IFf

        If_points = []

        for j in range(len(IF[0])):

            IF_temp = IF[:, j]
            If_points_temp = grid_points[IF_temp == 1]
            If_points.append(If_points_temp)

        return If_points

    def set_intrusion_frame_data(self, intrusion_frame_data, intrusion_network_points):

        """Adds the intrusion network points as coordinate 0 data for the intrusion frame

        Parameters
        ----------
        intrusion_frame_data : DataFrame,
                                intrusion frame data from model data
        intrusion_network_points: numpy array [x,y,z],
                                  outcome of self.create_intrusion_network()

        """
        # Coordinate 0 - Represents growth, isovalue 0 correspond to the intrusion network surface, gradient must be provided (ix,iy,iz):
        scaled_inet_points = intrusion_network_points[:, :3]
        coord_0_values = pd.DataFrame(scaled_inet_points, columns=["X", "Y", "Z"])
        coord_0_values["val"] = 0
        coord_0_values["coord"] = 0
        coord_0_values["feature_name"] = self.name
        coord_0_values["w"] = 1

        intrusion_frame_data_temp = pd.concat([intrusion_frame_data, coord_0_values])

        # if self.intrusion_network_gradients is not None:
        coord_0_data = intrusion_frame_data[intrusion_frame_data["coord"] == 0].copy()
        if len(coord_0_data) == 0:
            np.random.shuffle(self.intrusion_network_gradients)
            sliced_grads = self.intrusion_network_gradients[:200, :]

            coord_0_grads = pd.DataFrame(
                sliced_grads[:, :6], columns=["X", "Y", "Z", "gx", "gy", "gz"]
            )
            coord_0_grads["coord"] = 0
            coord_0_grads["feature_name"] = self.name
            coord_0_grads["w"] = 1

            intrusion_frame_data_complete = pd.concat(
                [intrusion_frame_data_temp, coord_0_grads]
            )

        else:
            intrusion_frame_data_complete = intrusion_frame_data_temp

        self.add_data_from_data_frame(intrusion_frame_data_complete)
        self.update_geometry(intrusion_frame_data_complete[["X", "Y", "Z"]].to_numpy())
        # return intrusion_frame_data

    def update(self):
        for i in range(3):
            self.builders[i].update()
