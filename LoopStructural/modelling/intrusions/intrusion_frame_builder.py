from ...modelling.features.builders import StructuralFrameBuilder
from ...modelling.features.fault import FaultSegment
from ...utils import getLogger, rng
from ...datatypes import BoundingBox

from typing import Union


logger = getLogger(__name__)

import numpy as np
import pandas as pd

try:
    from sklearn.cluster import KMeans
except ImportError as e:
    logger.error('Scikitlearn cannot be imported')
    raise e


class IntrusionFrameBuilder(StructuralFrameBuilder):
    def __init__(
        self,
        interpolatortype: Union[str, list],
        bounding_box: BoundingBox,
        nelements: Union[int, list] = 1000,
        model=None,
        **kwargs,
    ):
        """IntrusionBuilder set up the intrusion frame to build an intrusion
            The intrusion frame is curvilinear coordinate system of the intrusion, which is then used to compute the lateral and vertical exten of the intrusion.
            The object is constrained with the host rock geometry and other field measurements such as propagation and inflation vectors.

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        model : GeologicalModel
            reference to the model containing the fault
        """

        StructuralFrameBuilder.__init__(self, interpolatortype, bounding_box, nelements, **kwargs)

        self.origin = np.array([np.nan, np.nan, np.nan])
        self.maximum = np.array([np.nan, np.nan, np.nan])
        self.model = model
        self.minimum_origin = self.model.bounding_box[0, :]
        self.maximum_maximum = self.model.bounding_box[1, :]
        self.faults = []

        # -- intrusion frame parameters
        self.intrusion_network_contact = (
            None  # string, contact which is used to constrain coordinate 0 (roof/top or floor/base)
        )
        self.intrusion_other_contact = (
            None  # string, the contact which is NOT used to constrain coordinate 0
        )
        self.intrusion_network_type = (
            None  # string, (FAN: REMOVE shortest path option, no longer working)
        )
        self.intrusion_network_data = None
        self.other_contact_data = None
        self.grid_to_evaluate_ifx = np.zeros([1, 1])

        # -- host rock anisotropies exploited by the intrusion (lists of names and parameters)
        self.anisotropies_series_list = []
        self.anisotropies_series_parameters = {}
        self.anisotropies_fault_list = []
        self.anisotropies_fault_parameters = {}

        self.delta_contacts = None
        self.delta_faults = None
        self.intrusion_steps = None
        self.marginal_faults = None
        self.frame_c0_points = None
        self.frame_c0_gradients = None

        self.IFf = None
        self.IFc = None

    def update_geometry(self, points):
        self.origin = np.nanmin(np.array([np.min(points, axis=0), self.origin]), axis=0)
        self.maximum = np.nanmax(np.array([np.max(points, axis=0), self.maximum]), axis=0)
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

    def add_fault(self, fault: FaultSegment):
        """
        Add a fault to the geological feature builder

        Parameters
        ----------
        fault : FaultSegment
            A faultsegment to add to the geological feature

        Returns
        -------

        """

        self._up_to_date = False
        self.faults.append(fault)
        for i in range(3):
            self.builders[i].add_fault(fault)

    def set_intrusion_frame_c0_data(self, intrusion_data: pd.DataFrame):
        """
        Separate data between roof and floor contacts

        Parameters
        ----------
        Returns
        -------
        """
        other_contact = None
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

        self.intrusion_other_contact = other_contact
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

        if spacing is None:
            spacing = self.model.nsteps

        grid_points = self.model.regular_grid(nsteps=spacing, shuffle=False)

        self.grid_to_evaluate_ifx = grid_points

        return grid_points, spacing

    def add_contact_anisotropies(self, series_list: list = [], **kwargs):
        """
        Currently only used in 'Shortest path algorithm' (deprecated).
        Add to the intrusion network the anisotropies
        exploited by the intrusion (series-type geological features)
        Given a list of series-type features, this function evaluates contact points
        on each series and compute mean value and standard deviation.
        Different contacts of the same series are indentify using
        clustering algortihm. Mean and std deviation values will
        be used to identify each contact thoughout the model using
        the indicator functions.

        Parameters
        ----------
        series_list: list
            list of series-type features

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
            data_temp = self.intrusion_network_data[
                self.intrusion_network_data["intrusion_anisotropy"] == series.name
            ].copy()
            data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()
            series_i_vals = series.evaluate_value(data_array_temp)
            series_array = np.zeros((len(data_array_temp), 4))
            series_array[:, :3] = data_array_temp
            series_array[:, 3] = series_i_vals

            n_contacts = n_clusters[i]

            # -- use scalar field values to find different contacts
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

    def add_faults_anisotropies(self, fault_list: list = []):
        """
        Add to the intrusion network the anisotropies likely
        exploited by the intrusion (fault-type geological features)
        Given a list of fault features, evaluates the contact
        points on each fault and compute mean value and standard
        deviation. These values will be used to identify each
        fault with the indicator function.
        If no points in the fault, assign standard mean and std dev.

        Parameters
        ----------
        fault_list: list
            list of fault-type features

        Returns
        -------

        """
        if fault_list is not None:
            self.anisotropies_fault_list.append(fault_list)

            for i in range(len(fault_list)):
                if fault_list[i] in self.faults:  # remove pre-intrusion faults from faults list
                    self.faults.remove(fault_list[i])
                    for j in range(3):
                        self.builders[j].faults.remove(fault_list[i])

                fault = fault_list[i]
                data_temp = self.intrusion_network_data[
                    self.intrusion_network_data["intrusion_anisotropy"] == fault.name
                ].copy()
                data_array_temp = data_temp.loc[:, ["X", "Y", "Z"]].to_numpy()

                if data_temp.empty:
                    fault_i_mean = 0
                    fault_i_std = 0.1

                else:
                    # -- evaluate the value of the fault, evaluate_value called on a fault
                    # will return the scalar field value of the main structural feature
                    fault_i_vals = fault[0].evaluate_value(data_array_temp)
                    fault_i_mean = np.mean(fault_i_vals)
                    fault_i_std = np.std(fault_i_vals)

                self.anisotropies_fault_parameters[fault_list[i].name] = [
                    fault_list[i],
                    fault_i_mean,
                    fault_i_std,
                ]

    def set_intrusion_steps_parameters(self):
        """
        Use intrusion contact data (reference contact for intrusion frame) to compute the parameters related to each step.
        The parameters are the mean and std dev values of the stratigraphic units to both sides of the step.
        These parameters are then used to add more constraints to the intrusion framce coordinate 0

        For each step, at least one strigraphic unit and one fault must be provided.

        Parameters
        ----------

        Returns
        -------

        """

        if self.model.stratigraphic_column is None:
            logger.error(
                "Set stratigraphic column to model using model.set_stratigraphic_column(name of stratigraphic column)"
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
            unit_to_id = self.model.stratigraphic_column[series_to_name.name][unit_to_name].get(
                "id"
            )

            intrusion_network_data.loc[:, "model_values"] = self.model.evaluate_model(
                self.model.rescale(intrusion_network_data_xyz, inplace=False)
            )

            # -- check if step is within the same unit. If so, find clusters of data:
            if unit_from_name == unit_to_name:
                data_points_xyz = (
                    intrusion_network_data[intrusion_network_data["model_values"] == unit_from_id]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                series_values = series_from_name.evaluate_value(data_points_xyz)
                series_values_mod = series_values.reshape(len(series_values), 1)
                contact_clustering = KMeans(n_clusters=2, random_state=0).fit(series_values_mod)

                # contact 0
                z = np.ma.masked_not_equal(contact_clustering.labels_, 0)
                y = np.ma.masked_array(series_values, z.mask)
                contact_0_vals = np.ma.compressed(y)
                contact_0_mean = np.mean(contact_0_vals)
                contact_0_std = np.std(contact_0_vals)

                if contact_0_std == 0 or np.isnan(contact_0_std):
                    contact_0_std = std_backup

                # contact 1
                z = np.ma.masked_not_equal(contact_clustering.labels_, 1)
                y = np.ma.masked_array(series_values, z.mask)
                contact_1_vals = np.ma.compressed(y)
                contact_1_mean = np.mean(contact_1_vals)
                contact_1_std = np.std(contact_1_vals)

                if contact_1_std == 0 or np.isnan(contact_1_std):
                    contact_1_std = std_backup

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

            else:  # -- step between different stratigraphic units
                data_points_from_xyz = (
                    intrusion_network_data[intrusion_network_data["model_values"] == unit_from_id]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                step_structure_points_vals = step_structure[0].evaluate_value(data_points_from_xyz)
                if len(data_points_from_xyz) == 0:  # no data points in strat unit
                    unit_from_min = self.model.stratigraphic_column[series_from_name.name][
                        unit_from_name
                    ].get("min")
                    unit_from_max = self.model.stratigraphic_column[series_from_name.name][
                        unit_from_name
                    ].get("max")
                    self.intrusion_steps[step_i]["unit_from_mean"] = (
                        unit_from_min + (unit_from_max - unit_from_min) / 2
                    )
                    self.intrusion_steps[step_i]["unit_from_std"] = std_backup
                else:
                    series_values = series_from_name.evaluate_value(data_points_from_xyz)
                    mask = step_structure_points_vals < 0
                    if len(mask) > 0:
                        series_values_mod = np.ma.compressed(
                            np.ma.masked_array(series_values, step_structure_points_vals < 0)
                        )
                    else:
                        series_values_mod = series_values
                    step["unit_from_mean"] = np.nanmean(series_values_mod)
                    step["unit_from_std"] = np.nanstd(series_values_mod)

                if step["unit_from_std"] == 0:
                    step["unit_from_std"] = std_backup

                data_points_to_xyz = (
                    intrusion_network_data[intrusion_network_data["model_values"] == unit_to_id]
                    .loc[:, ["X", "Y", "Z"]]
                    .copy()
                    .to_numpy()
                )
                step_structure_points_vals = step_structure[0].evaluate_value(data_points_to_xyz)
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
                    mask = step_structure_points_vals > 0
                    if len(mask) > 0:
                        series_values_mod = np.ma.compressed(
                            np.ma.masked_array(series_values, step_structure_points_vals > 0)
                        )
                    else:
                        series_values_mod = series_values
                    step["unit_to_mean"] = np.nanmean(series_values_mod)
                    step["unit_to_std"] = np.nanstd(series_values_mod)
                    check_mean = step["unit_to_mean"]
                    check_std = step["unit_to_std"]

                    if np.isnan(check_mean):
                        step["unit_to_mean"] = 40

                    if np.isnan(check_std):
                        step["unit_to_std"] = std_backup

                if step["unit_to_std"] == 0:
                    step["unit_to_std"] = std_backup

    def set_marginal_faults_parameters(self):
        """
        Use intrusion contact data (reference contact for intrusion frame) to compute the parameters related to each fault.
        The parameters are the mean and std dev values of the stratigraphic units to the hanging wall or foot wall of the fault.
        These parameters are then used to add more constraints to the intrusion framce coordinate 0

        For each fault, at least one strigraphic unit and one fault must be provided.

        Parameters
        ----------

        Returns
        -------

        """

        # set data
        intrusion_frame_c0_data = self.intrusion_network_data.copy()
        intrusion_frame_c0_data_xyz = (
            intrusion_frame_c0_data.loc[:, ["X", "Y", "Z"]].copy().to_numpy()
        )
        std_backup = 25

        for fault_i in self.marginal_faults.keys():
            marginal_fault = self.marginal_faults[fault_i].get("structure")
            block = self.marginal_faults[fault_i].get("block")  # hanging wall or foot wall
            self.marginal_faults[fault_i].get("emplacement_mechanism")
            series_name = self.marginal_faults[fault_i].get("series")

            series_values_temp = series_name.evaluate_value(intrusion_frame_c0_data_xyz)
            faults_values_temp = marginal_fault[0].evaluate_value(intrusion_frame_c0_data_xyz)

            if block == "hanging wall":
                series_values = series_values_temp[faults_values_temp > 0]

            elif block == "foot wall":
                series_values = series_values_temp[faults_values_temp < 0]

            self.marginal_faults[fault_i]["series_mean"] = np.mean(series_values)

            series_std = np.std(series_values)

            if series_std == 0:
                series_std = std_backup

            self.marginal_faults[fault_i]["series_std"] = series_std
            self.marginal_faults[fault_i]["series_vals"] = series_values

    def set_intrusion_frame_parameters(
        self, intrusion_data: pd.DataFrame, intrusion_frame_parameters: dict, **kwargs
    ):
        """
        Set variables to create intrusion network.

        Parameters
        ----------
        intrusion_data: pd.DataFrame
            intrusion contact data
        intrusion_network_input: dict
            contact : string, contact of the intrusion to be used to create the network (roof or floor)
            type : string, type of algorithm to create the intrusion network (interpolated or shortest path).
                    Shortest path is recommended when intrusion contact is not well constrained
            contacts_anisotropies : list of series-type features involved in intrusion emplacement
            structures_anisotropies : list of fault-type features involved in intrusion emplacement
            sequence_anisotropies : list of anisotropies to look for the shortest path. It could be only starting and end point.

        Returns
        -------

        """

        self.intrusion_network_contact = intrusion_frame_parameters.get("contact", "floor")

        self.set_intrusion_frame_c0_data(intrusion_data)  # separates roof and floor data

        # if self.intrusion_network_type == "interpolated":

        self.gradients_constraints_weight = intrusion_frame_parameters.get("g_w", None)

        intrusion_steps = intrusion_frame_parameters.get("intrusion_steps", None)

        marginal_faults = intrusion_frame_parameters.get("marginal_faults", None)

        if intrusion_steps is not None:
            self.intrusion_steps = intrusion_steps

            self.delta_contacts = intrusion_frame_parameters.get(
                "delta_c", [1] * len(intrusion_steps)
            )

            self.delta_faults = intrusion_frame_parameters.get(
                "delta_f", [1] * len(intrusion_steps)
            )

            self.set_intrusion_steps_parameters()  # function to compute steps parameters

            fault_anisotropies = []
            for step in self.intrusion_steps.keys():
                fault_anisotropies.append(self.intrusion_steps[step].get("structure"))

            self.add_faults_anisotropies(fault_anisotropies)

        if marginal_faults is not None:
            self.marginal_faults = marginal_faults

            self.delta_contacts = intrusion_frame_parameters.get(
                "delta_c", [1] * len(marginal_faults)
            )

            self.delta_faults = intrusion_frame_parameters.get(
                "delta_f", [1] * len(marginal_faults)
            )

            self.set_marginal_faults_parameters()

            fault_anisotropies = []
            for fault in self.marginal_faults.keys():
                fault_anisotropies.append(self.marginal_faults[fault].get("structure"))

                self.add_faults_anisotropies(fault_anisotropies)

        # add contact anisotropies list
        contact_anisotropies = intrusion_frame_parameters.get("contact_anisotropies", None)

        self.anisotropies_series_list = contact_anisotropies

    def indicator_function_contacts(self, delta=None):
        """
        Function only used for Shortest Path method (Deprecated)
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

            # series_id.faults_enabled = True
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

    def create_constraints_for_c0(self, **kwargs):
        """
        Creates a numpy array containing (x,y,z) coordinates of synthetic points located
        in the intrusion reference contact.
        The intrusion reference contact is the one used to build c0 of the intrusion frame.
        --- Currently only works if steps OR marginal faults are present ---

        Parameters
        ----------

        Returns
        -------
        intrusion_contact_points = numpy array
        """

        # --- data points of intrusion reference contact (roof or floor)
        inet_points_xyz = self.intrusion_network_data.loc[:, ["X", "Y", "Z"]].to_numpy()

        intrusion_reference_contact_points = inet_points_xyz

        grid_points, spacing = self.create_grid_for_indicator_fxs()

        # --- more constraints if steps or marginal fault is present:
        if self.intrusion_steps is not None:
            If = self.indicator_function_faults(delta=self.delta_faults)

            if len(np.where(If == 1)[0]) == 0:
                logger.error("No faults identified, you may increase the value of delta_f")

            If_sum = np.sum(If, axis=1)

            # -- evaluate grid points in series

            for i, step_i in enumerate(self.intrusion_steps.keys()):
                delta_contact = self.delta_contacts[i]
                if isinstance(delta_contact, list):
                    delta_contact0 = delta_contact[0]
                    delta_contact1 = delta_contact[1]
                else:
                    delta_contact0 = delta_contact
                    delta_contact1 = delta_contact

                step_fault = self.intrusion_steps[step_i].get("structure")
                series_from_name = self.intrusion_steps[step_i].get("series_from")
                series_to_name = self.intrusion_steps[step_i].get("series_to")

                fault_gridpoints_vals = step_fault[0].evaluate_value(grid_points)

                series_from_gridpoints_vals = series_from_name.evaluate_value(grid_points)

                if series_from_name == series_to_name:
                    series_to_gridpoints_vals = series_from_gridpoints_vals

                else:
                    series_to_gridpoints_vals = series_to_name.evaluate_value(grid_points)

                contacts0_val_min = self.intrusion_steps[step_i].get("unit_from_mean") - (
                    self.intrusion_steps[step_i].get("unit_from_std") * delta_contact0
                )
                contacts0_val_max = self.intrusion_steps[step_i].get("unit_from_mean") + (
                    self.intrusion_steps[step_i].get("unit_from_std") * delta_contact0
                )
                contacts1_val_min = self.intrusion_steps[step_i].get("unit_to_mean") - (
                    self.intrusion_steps[step_i].get("unit_to_std") * delta_contact1
                )
                contacts1_val_max = self.intrusion_steps[step_i].get("unit_to_mean") + (
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

                step_i_constraints_temp = grid_points[
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

                region = self.intrusion_steps[step_i].get("region", None)
                if region is None:
                    step_i_constraints = step_i_constraints_temp
                else:
                    mask = region(step_i_constraints_temp)
                    step_i_constraints = step_i_constraints_temp[mask]

                intrusion_reference_contact_points = np.vstack(
                    [intrusion_reference_contact_points, step_i_constraints]
                )

                splits_from_sill_name = self.intrusion_steps[step_i].get("splits_from", None)

                # check if sill comes from anothe sill
                # (add all the constraint from original sill, this have to be changed and adapted so it adds constraints for specific faults)

                if splits_from_sill_name is not None:
                    splits_from_sill_steps = self.model.__getitem__(
                        splits_from_sill_name
                    ).intrusion_frame.builder.intrusion_steps
                    for step_j in splits_from_sill_steps.keys():
                        step_j_hg_constraints = splits_from_sill_steps[step_j].get("constraints_hw")
                        intrusion_reference_contact_points = np.vstack(
                            [intrusion_reference_contact_points, step_j_hg_constraints]
                        )

            # --- more constraints if steps or marginal fault is present:
        if self.marginal_faults is not None:
            If = self.indicator_function_faults(delta=self.delta_faults)

            if len(np.where(If == 1)[0]) == 0:
                logger.error("No faults identified, yo may increase value of delta_f")

            If_sum = np.sum(If, axis=1)

            # evaluate grid points in series
            for fault_i in self.marginal_faults.keys():
                delta_contact = self.marginal_faults[fault_i].get("delta_c", 1)
                marginal_fault = self.marginal_faults[fault_i].get("structure")
                block = self.marginal_faults[fault_i].get("block")  # hanging wall or foot wall
                self.marginal_faults[fault_i].get("emplacement_mechanism")
                series_name = self.marginal_faults[fault_i].get("series")

                fault_gridpoints_vals = marginal_fault[0].evaluate_value(grid_points)
                series_gridpoints_vals = series_name.evaluate_value(grid_points)

                contact_min = self.marginal_faults[fault_i].get("series_mean") - (
                    self.marginal_faults[fault_i].get("series_std") * delta_contact
                )
                contact_max = self.marginal_faults[fault_i].get("series_mean") + (
                    self.marginal_faults[fault_i].get("series_std") * delta_contact
                )

                if block == "hanging wall":
                    marginalfault_i_constraints_temp = grid_points[
                        np.logical_and(
                            If_sum == 0,
                            np.logical_and(
                                fault_gridpoints_vals >= 0,
                                np.logical_and(
                                    series_gridpoints_vals >= contact_min,
                                    series_gridpoints_vals <= contact_max,
                                ),
                            ),
                        )
                    ]

                elif block == "foot wall":
                    marginalfault_i_constraints_temp = grid_points[
                        np.logical_and(
                            If_sum == 0,
                            np.logical_and(
                                fault_gridpoints_vals <= 0,
                                np.logical_and(
                                    series_gridpoints_vals >= contact_min,
                                    series_gridpoints_vals <= contact_max,
                                ),
                            ),
                        )
                    ]

                region = self.marginal_faults[fault_i].get("region", None)

                if region is None:
                    marginalfault_i_constraints = marginalfault_i_constraints_temp
                else:
                    mask = region(marginalfault_i_constraints_temp)
                    marginalfault_i_constraints = marginalfault_i_constraints_temp[mask]

                intrusion_reference_contact_points = np.vstack(
                    [intrusion_reference_contact_points, marginalfault_i_constraints]
                )

        self.frame_c0_points = intrusion_reference_contact_points

        # -- Gradient contraints for intrusion frame c0
        # currently only used if no inflation data, and for one/main series anisotropy
        # Evaluate points in gradient of stratigraphy, and exclude points around faults

        series_id = self.anisotropies_series_list[0]
        stratigraphy_gradient_grid_points = series_id.evaluate_gradient(grid_points)

        # If intrusion frame c0 is built usinf roof/top contact, then change vector direction
        if self.intrusion_network_contact == "roof" or self.intrusion_network_contact == "top":
            grid_points_inflation = stratigraphy_gradient_grid_points * (-1)
        else:
            grid_points_inflation = stratigraphy_gradient_grid_points

        grid_points_and_inflation_all = np.hstack([grid_points, grid_points_inflation])

        if self.intrusion_steps is None:
            If_sum = np.zeros(
                len(grid_points_and_inflation_all)
            ).T  # mask for region without faults afecting the intrusion

        grid_points_and_inflation = grid_points_and_inflation_all[If_sum == 0]
        self.frame_c0_gradients = grid_points_and_inflation

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

    def set_intrusion_frame_data(self, intrusion_frame_data):  # , intrusion_network_points):
        """Adds the intrusion network points as coordinate 0 data for the intrusion frame

        Parameters
        ----------
        intrusion_frame_data : DataFrame,
                                intrusion frame data from model data

        """
        # Coordinate 0 - Represents growth, isovalue 0 represents the location of the roof or floor contact:
        c0_points = self.frame_c0_points[:, :3]
        coord_0_values = pd.DataFrame(c0_points, columns=["X", "Y", "Z"])
        coord_0_values["val"] = 0
        coord_0_values["coord"] = 0
        coord_0_values["feature_name"] = self.name
        coord_0_values["w"] = 1

        intrusion_frame_data_temp = pd.concat([intrusion_frame_data, coord_0_values])

        coord_0_data = intrusion_frame_data[intrusion_frame_data["coord"] == 0].copy()
        if len(coord_0_data) == 0:
            rng.shuffle(self.frame_c0_gradients)
            if self.gradients_constraints_weight is None:
                n_grad_constraints = 100
            else:
                n_grad_constraints = int(
                    len(self.frame_c0_gradients) * self.gradients_constraints_weight
                )
            sliced_grads = self.frame_c0_gradients[:n_grad_constraints, :]

            coord_0_grads = pd.DataFrame(
                sliced_grads[:, :6], columns=["X", "Y", "Z", "gx", "gy", "gz"]
            )
            coord_0_grads["coord"] = 0
            coord_0_grads["feature_name"] = self.name
            coord_0_grads["w"] = 1

            intrusion_frame_data_complete = pd.concat([intrusion_frame_data_temp, coord_0_grads])

        else:
            intrusion_frame_data_complete = intrusion_frame_data_temp

        self.add_data_from_data_frame(intrusion_frame_data_complete)
        self.update_geometry(intrusion_frame_data_complete[["X", "Y", "Z"]].to_numpy())

    def update(self):
        for i in range(3):
            self.builders[i].update()
