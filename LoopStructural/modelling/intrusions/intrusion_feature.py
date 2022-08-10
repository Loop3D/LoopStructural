import numpy as np
import pandas as pd

# import logging
from LoopStructural.utils import getLogger
from LoopStructural.modelling.features import FeatureType

logger = getLogger(__name__)


class IntrusionFeature:
    """
    Intrusion feature is a class to represent an intrusion, using a distance scalar field to the intrusion contact.
    Threshold distances are simulated along the intrusion frame coordinates, and simulation are constrained
    with conceptual geometrical model of the lateral and vertical intrusion extent.

    """

    def __init__(self, frame, builder, name="UnnamedIntrusion", model=None):

        """
        Parameters
        ----------
        name: string

        model: GeologicalModel

        Returns
        ----------
        intrusion_feature :  IntrusionFeature

        """

        self.name = name
        self.model = model
        self.intrusion_frame = frame
        self.builder = builder
        self.type = FeatureType.INTRUSION
        # simulated thresholds:
        self._lateral_simulated_thresholds = None
        self._growth_simulated_thresholds = None

    @property
    def lateral_simulated_thresholds(self):
        self.builder.up_to_date()

        return self._lateral_simulated_thresholds

    @lateral_simulated_thresholds.setter
    def lateral_simulated_thresholds(self, lateral_simulated_thresholds):
        # TODO check type is correct and will work?
        self._lateral_simulated_thresholds = lateral_simulated_thresholds

    @property
    def growth_simulated_thresholds(self):
        self.builder.up_to_date()

        return self._growth_simulated_thresholds

    @growth_simulated_thresholds.setter
    def growth_simulated_thresholds(self, growth_simulated_threshold):
        # TODO check type is correct and will work?
        self._growth_simulated_thresholds = growth_simulated_threshold

    @property
    def lateral_sgs_input_data(self):
        self.builder.up_to_date()
        return self.builder.lateral_sgs_input_data

    @property
    def vertical_sgs_input_data(self):
        self.builder.up_to_date()
        return self.builder.vertical_sgs_input_data

    def min(self):

        if self.model is None:
            return 0
        return np.nanmin(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def max(self):
        """
        Calculate average of the support values

        Returns
        -------
        max : double
            max value of the feature evaluated on a (10,10,10) grid in model area
        """
        if self.model is None:
            return 0
        return np.nanmax(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def set_intrusion_builder(self, builder):
        self.intrusion_builder = builder

    def set_intrusion_frame(self, intrusion_frame):
        self.intrusion_feature_frame = intrusion_frame

    def set_model(self, model):
        self.model = model

    def evaluate_value(self, points):

        """
        Computes a distance scalar field to the intrusion contact (isovalue = 0).

        Parameters
        ------------
        points : numpy array (x,y,z),  points where the IntrusionFeature is evaluated.

        Returns
        ------------
        intrusion_sf : numpy array, contains distance to intrusion contact

        """
        self.builder.up_to_date()

        simulation_g_data = self.growth_simulated_thresholds
        simulation_s_data = self.lateral_simulated_thresholds

        # ---> returns indicator function and scalar field with isovalue 0 = intrusion boundary

        # compute coordinates values for each evaluated point
        intrusion_coord0_pts = self.intrusion_frame[0].evaluate_value(points)
        intrusion_coord1_pts = self.intrusion_frame[1].evaluate_value(points)
        intrusion_coord2_pts = self.intrusion_frame[2].evaluate_value(points)

        # ------ lateral extent thresholds for each of the evaluated points -------------

        # if lateral extent values were simulated
        if simulation_s_data is None:
            print("No simulation for lateral extent")
        else:
            # print("Assigning lateral thresholds")
            simulation_s_data.sort_values(["coord1"], ascending=[True], inplace=True)

            # containers for thresholds
            s_minside_threshold = np.zeros(len(intrusion_coord1_pts))
            s_maxside_threshold = np.zeros(len(intrusion_coord1_pts))

            # simulated values (datframe to array)
            simulation_s_data_coord1 = simulation_s_data["coord1"].to_numpy()
            simulated_smin_values = simulation_s_data["min_l_threshold"].to_numpy()
            simulated_smax_values = simulation_s_data["max_l_threshold"].to_numpy()

            # find index of closest value to each point being evaluated, and assign simulated s thresholds
            pi_p = abs(
                intrusion_coord1_pts[:, None] - simulation_s_data_coord1[None, :]
            )
            indexS = np.argmin(pi_p, axis=1)

            for i in range(len(intrusion_coord0_pts)):
                s_minside_threshold[i] = simulated_smin_values[indexS[i]]
                s_maxside_threshold[i] = simulated_smax_values[indexS[i]]

            # ---- indicator function to define volume of intrusion - w/out considering growth threshold
            indicator_fxS_boolean = (s_minside_threshold <= intrusion_coord2_pts) * (
                intrusion_coord2_pts <= s_maxside_threshold
            )
            indicator_fxS = indicator_fxS_boolean.astype("int64")

        if simulation_g_data is None:
            print("No simulation for vertical extent")
        else:
            # containers for thresholds
            # print("Assigning vertical thresholds")
            g_minside_threshold = np.zeros(len(intrusion_coord1_pts))
            g_maxside_threshold = np.zeros(len(intrusion_coord1_pts))

            # simulated values (dataframe to array)
            simulation_g_data_coord1 = simulation_g_data["coord1"].to_numpy()
            simulation_g_data_coord2 = simulation_g_data["coord2"].to_numpy()
            simulated_gmin_values = simulation_g_data["g_minimum"].to_numpy()
            simulated_gmax_values = simulation_g_data["g_maximum"].to_numpy()

            # find index of closest value to each point being evaluated, and assign simulated s thresholds
            pi_p = abs(
                intrusion_coord1_pts[:, None] - simulation_g_data_coord1[None, :]
            )  # p_points - p_simulated_data
            si_s = abs(
                intrusion_coord2_pts[:, None] - simulation_g_data_coord2[None, :]
            )

            s_min = np.around(
                si_s / (si_s.min(axis=1)[:, None] + np.finfo("float").eps), 2
            )  # flag with 1 the minimum value
            p_min = np.around(
                pi_p / (pi_p.min(axis=1)[:, None] + np.finfo("float").eps), 2
            )  # flag with 1 the minimum value

            indexG = np.argmin(abs(1 - (s_min * p_min)), axis=1)

            for i in range(len(intrusion_coord0_pts)):
                g_minside_threshold[i] = simulated_gmin_values[indexG[i]]
                g_maxside_threshold[i] = simulated_gmax_values[indexG[i]]

            indicator_fxG_boolean = (g_minside_threshold <= intrusion_coord0_pts) * (
                intrusion_coord0_pts <= g_maxside_threshold
            )
            indicator_fxG = indicator_fxG_boolean.astype("int64")

        indicator_fx_boolean = indicator_fxS_boolean * indicator_fxG_boolean
        indicator_fx = indicator_fx_boolean.astype("int64")

        #         ------- intrusion_sf: final distance scalar field
        # Transform the scalar fields given by the frame coordinates, using the thresholds.
        # This aims to generate a scalar field with its isovalue = 0 on the intrusion contact

        mid_point = g_minside_threshold + (
            (g_maxside_threshold - g_minside_threshold) / 2
        )

        a = intrusion_coord2_pts >= s_maxside_threshold
        b = intrusion_coord2_pts <= s_minside_threshold
        c = (
            (s_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < s_maxside_threshold)
            * (intrusion_coord0_pts <= g_minside_threshold)
        )
        d = (
            (s_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < s_maxside_threshold)
            * (intrusion_coord0_pts >= g_maxside_threshold)
        )
        e = (
            (s_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < s_maxside_threshold)
            * (mid_point >= intrusion_coord0_pts)
            * (intrusion_coord0_pts > g_minside_threshold)
        )
        f = (
            (s_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < s_maxside_threshold)
            * (mid_point < intrusion_coord0_pts)
            * (intrusion_coord0_pts < g_maxside_threshold)
        )

        mod_Smin_thresholds = intrusion_coord2_pts - s_minside_threshold
        mod_Smax_thresholds = intrusion_coord2_pts - s_maxside_threshold
        mod_Gmin_thresholds = intrusion_coord0_pts - g_minside_threshold
        mod_Gmax_thresholds = intrusion_coord0_pts - g_maxside_threshold

        intrusion_sf = (
            a * mod_Smax_thresholds
            + b * abs(mod_Smin_thresholds)
            + c * abs(mod_Gmin_thresholds)
            + d * mod_Gmax_thresholds
            - e * mod_Gmin_thresholds
            + f * mod_Gmax_thresholds
        )

        self.evaluated_points = [
            points,
            intrusion_coord0_pts,
            intrusion_coord1_pts,
            intrusion_coord2_pts,
        ]
        self.intrusion_indicator_function = indicator_fx

        return intrusion_sf
