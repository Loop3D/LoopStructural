import numpy as np
import pandas as pd

# from LoopStructural import GeologicalModel

import logging

import numpy as np

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class IntrusionFeature:
    """
    Intrusion feature is a class to represent an intrusion, using a distance scalar field to the intrusion contact.
    """

    def __init__(self, name, structural_frame=None, model=None):

        """
        Parameters
        ----------
        name: string
        simulation_gdata: dataframe containing thresholds distances to constraint lateral extent of intrusion
        simulation_sdata: dataframe containing thresholds distances to constraint vertical extent of intrusion
        inet_plygon: only if lateral extent is provided using a polygon
        structural_frame: StructuralFrame
        model: GeologicalModel
        """

        self.name = name
        self.structural_frame = structural_frame
        self.model = model

        self.simulation_growth_data = None
        self.simulation_lateral_data = None
        # self.inet_poly = inet_polygon

        self.intrusion_feature_network = None
        self.intrusion_feature_frame = None
        self.intrusion_feature_body = None

        self.intrusion_indicator_function = None
        self.evaluated_points = None

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

    def set_simulation_growth_data(self, simulation_gdata):
        self.simulation_growth_data = simulation_gdata

    def set_simulation_lateral_data(self, simulation_sdata):
        self.simulation_lateral_data = simulation_sdata

    def set_intrusion_network(self, intrusion_network):
        self.intrusion_feature_network = intrusion_network

    def set_intrusion_frame(self, intrusion_frame):
        self.intrusion_feature_frame = intrusion_frame

    def set_intrusion_body(self, intrusion_body):
        self.intrusion_feature_body = intrusion_body

    def evaluate_value(self, points):

        """
        points : array (x,y,z),  points where the intrusion is evaluated.

        """

        simulation_g_data = self.simulation_growth_data
        simulation_s_data = self.simulation_lateral_data
        inet_polygon = None
        intrusion_frame = self.structural_frame
        model = self.model

        # ---> returns indicator function and scalar field with isovalue 0 = intrusion boundary

        # compute coordinates values for each evaluated point
        intrusion_coord0_pts = intrusion_frame[0].evaluate_value(points)
        intrusion_coord1_pts = intrusion_frame[1].evaluate_value(points)
        intrusion_coord2_pts = intrusion_frame[2].evaluate_value(points)

        # ------ lateral extent thresholds for each of the evaluated points -------------

        # if lateral extent values were simulated
        if simulation_s_data is None:
            print("No simultion for lateral extent")
        else:
            print("asssigning lateral thresholds")
            simulation_s_data.sort_values(["coord1"], ascending=[True], inplace=True)

            # containers for thresholds
            s_minside_threshold = np.zeros(len(intrusion_coord1_pts))
            s_maxside_threshold = np.zeros(len(intrusion_coord1_pts))

            # simulated values (datframe to array)
            simulation_s_data_coord1 = simulation_s_data["coord1"].to_numpy()
            simulated_smin_values = simulation_s_data["min_s_threshold"].to_numpy()
            simulated_smax_values = simulation_s_data["max_s_threshold"].to_numpy()

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
            print("No simultion for vertical extent")
        else:
            # containers for thresholds
            print("asssigning vertical thresholds")
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
                si_s / si_s.min(axis=1)[:, None], 2
            )  # flag with 1 the minimum value
            p_min = np.around(
                pi_p / pi_p.min(axis=1)[:, None], 2
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
        # print('modifying thresholds')

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
            * (0 >= intrusion_coord0_pts)
            * (intrusion_coord0_pts > g_minside_threshold)
        )
        f = (
            (s_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < s_maxside_threshold)
            * (0 < intrusion_coord0_pts)
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
    
        return intrusion_sf  # , indicator_fx, intrusion_sf2_temp, intrusion_sf0_temp,

    def evaluate_value2(self, points):

        """
        points : array (x,y,z),  points where the intrusion is evaluated.

        """

        simulation_g_data = self.simulation_growth_data
        simulation_s_data = self.simulation_lateral_data
        inet_polygon = None
        intrusion_frame = self.structural_frame
        model = self.model

        # points = model.scale(points, inplace = False)

        # ---> returns indicator function and scalar field with isovalue 0 = intrusion boundary

        # compute coordinates values for each evaluated point
        # intrusion_coord0_pts = intrusion_frame[0].evaluate_value(model.rescale(points, inplace = False))
        # intrusion_coord1_pts = intrusion_frame[1].evaluate_value(model.rescale(points, inplace = False))
        # intrusion_coord2_pts = intrusion_frame[2].evaluate_value(model.rescale(points, inplace = False))
        intrusion_coord0_pts = intrusion_frame[0].evaluate_value(points)
        intrusion_coord1_pts = intrusion_frame[1].evaluate_value(points)
        intrusion_coord2_pts = intrusion_frame[2].evaluate_value(points)
        # ------ lateral extent thresholds for each of the evaluated points -------------
        indicator_fx_intrusion = np.zeros([len(points)])

        # if lateral extent values were simulated
        if simulation_s_data is not None:

            simulation_s_data.sort_values(["coord1"], ascending=[True], inplace=True)
            s_minside_threshold = np.zeros(len(intrusion_coord1_pts))
            s_maxside_threshold = np.zeros(len(intrusion_coord1_pts))

            for i in range(len(simulation_s_data)):
                for j in range(len(intrusion_coord0_pts)):
                    if s_minside_threshold[j] == 0 and s_maxside_threshold[j] == 0:
                        if (
                            intrusion_coord1_pts[j]
                            <= simulation_s_data.loc[i, "coord1"]
                        ):
                            s_minside_threshold[j] = simulation_s_data.loc[
                                i, "min_s_threshold"
                            ]
                            s_maxside_threshold[j] = simulation_s_data.loc[
                                i, "max_s_threshold"
                            ]

            # ---- indicator function to define area of intrusion - w/out considering growth threshold
            min_inet_coord0 = simulation_g_data.loc[:, "g_minimum"].min()
            for i in range(len(points)):
                if indicator_fx_intrusion[i] == 0:
                    if (
                        intrusion_coord0_pts[i] >= min_inet_coord0
                    ):  # >=0 #constaint for intrusion network
                        if (
                            s_minside_threshold[i]
                            <= intrusion_coord2_pts[i]
                            <= s_maxside_threshold[i]
                        ):  # constraint for lateral extent
                            indicator_fx_intrusion[i] = 1

        # if lateral extent polygon is given
        elif inet_polygon != None:
            from shapely.geometry import Point, Polygon

            for i in range(len(intrusion_coord1_pts)):
                p = Point(intrusion_coord1_pts[i], intrusion_coord2_pts[i])
                if (
                    intrusion_coord0_pts[i] >= min_inet_coord0
                ):  # >=0 #constaint for intrusion network
                    if p.within(inet_polygon) == True:
                        indicator_fx_intrusion[i] = 1

        # ------ growth threshold (only computed for points indicator_fx_intrusion[i] == 1)
        growth_threshold = np.zeros([len(intrusion_coord0_pts)])
        inet_threshold = np.zeros([len(intrusion_coord0_pts)])
        for i in range(len(points)):
            if indicator_fx_intrusion[i] == 1:
                for j in range(len(simulation_g_data)):
                    if intrusion_coord1_pts[i] <= simulation_g_data.loc[j, "coord1"]:
                        if (
                            intrusion_coord2_pts[i]
                            <= simulation_g_data.loc[j, "coord2"]
                        ):
                            growth_threshold[i] = simulation_g_data.loc[j, "g_maximum"]
                            inet_threshold[i] = simulation_g_data.loc[j, "g_minimum"]
                            break

        # ------ Correct indicator function using growth threshold

        for i in range(len(points)):
            if indicator_fx_intrusion[i] == 1:
                if (
                    intrusion_coord0_pts[i] > growth_threshold[i]
                    or intrusion_coord0_pts[i] < inet_threshold[i]
                ):
                    indicator_fx_intrusion[i] = 0

        #         ------- intrusion_sf: final distance scalar field

        # containers for distance scalar field of each coordinate
        intrusion_sf0_temp = np.zeros(len(intrusion_coord0_pts))
        intrusion_sf1_temp = np.zeros(len(intrusion_coord1_pts))
        intrusion_sf2_temp = np.zeros(len(intrusion_coord2_pts))
        intrusion_sf = np.zeros(len(points))

        for i in range(len(points)):  # growth scalar field
            if intrusion_coord0_pts[i] < 0:
                intrusion_sf0_temp[i] = abs(intrusion_coord0_pts[i])
            elif intrusion_coord0_pts[i] >= growth_threshold[i]:
                intrusion_sf0_temp[i] = abs(
                    intrusion_coord0_pts[i] - growth_threshold[i]
                )
            else:
                intrusion_sf0_temp[i] = intrusion_coord0_pts[i] - growth_threshold[i]

        # for i in range(len(points)): #strike scalar field
        #     if intrusion_coord2_pts[i] >= s_maxside_threshold[i]:
        #         intrusion_sf2_temp[i] = intrusion_coord2_pts[i] - s_maxside_threshold[i]
        #     elif intrusion_coord2_pts[i] <= s_minside_threshold[i]:
        #         intrusion_sf2_temp[i] = abs(intrusion_coord2_pts[i] - s_minside_threshold[i])
        #     elif 0 <= intrusion_coord2_pts[i] < s_maxside_threshold[i]:
        #         intrusion_sf2_temp[i] = intrusion_coord2_pts[i] - s_maxside_threshold[i]
        #     elif 0 > intrusion_coord2_pts[i] > s_minside_threshold[i]:
        #         intrusion_sf2_temp[i] = -(intrusion_coord2_pts[i] - s_minside_threshold[i])
        #     else: continue

        intrusion_sf = intrusion_sf0_temp

        # for i in range(len(indicator_fx_intrusion)):
        #     if indicator_fx_intrusion[i] == 1:
        #         intrusion_sf[i] = -10
        #     else:
        #         intrusion_sf[i] = 1 + abs(intrusion_coord0_pts[i]) + abs(intrusion_coord1_pts[i]) + abs(intrusion_coord2_pts[i])

        self.evaluated_points = [
            points,
            intrusion_coord0_pts,
            intrusion_coord1_pts,
            intrusion_coord2_pts,
        ]
        self.intrusion_indicator_function = indicator_fx_intrusion
        return intrusion_sf  # , indicator_fx_intrusion
