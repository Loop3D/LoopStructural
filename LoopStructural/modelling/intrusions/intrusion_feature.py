from typing import Optional
import numpy as np
import pandas as pd
from LoopStructural.modelling.features import BaseFeature
from LoopStructural.modelling.features import FeatureType

# import logging
from ...utils import getLogger
from scipy.interpolate import Rbf

logger = getLogger(__name__)


class IntrusionFeature(BaseFeature):
    """
    Intrusion feature is a class to represent an intrusion, using a distance scalar field to the intrusion contact.
    Threshold distances are simulated along the intrusion frame coordinates, and simulation are constrained
    with conceptual geometrical model of the lateral and vertical intrusion extent.

    """

    def __init__(
        self,
        frame,
        builder,
        # faults=[],
        name="UnnamedIntrusion",
        model=None,
    ):
        """
        Parameters
        ----------
        name: string

        model: GeologicalModel

        Returns
        ----------
        intrusion_feature :  IntrusionFeature

        """

        BaseFeature.__init__(self, name=name, builder=builder, model=model)

        self.intrusion_frame = frame
        self.type = FeatureType.INTRUSION
        self.assisting_faults = {}

    def set_intrusion_frame(self, intrusion_frame):
        self.intrusion_feature_frame = intrusion_frame

    def set_model(self, model):
        self.model = model

    def add_assisting_faults(self, faults_dictionary):

        self.assisting_faults = faults_dictionary

    def interpolate_lateral_thresholds(self, points_coord1):

        conceptual_model = self.builder.lateral_extent_model
        inputsimdata_minL = self.builder.data_for_lateral_extent_calculation[0]
        inputsimdata_maxL = self.builder.data_for_lateral_extent_calculation[1]

        minL_inputdata_coord1 = inputsimdata_minL.coord1.to_numpy()
        # minL_inputdata_coord2 = inputsimdata_minL.coord2.to_numpy()
        minL_inputdata_residual = inputsimdata_minL.l_residual.to_numpy()
        # minL_inputdata_conceptual = inputsimdata_minL.l_conceptual.to_numpy()

        maxL_inputdata_coord1 = inputsimdata_maxL.coord1.to_numpy()
        # maxL_inputdata_coord2 = inputsimdata_maxL.coord2.to_numpy()
        maxL_inputdata_residual = inputsimdata_maxL.l_residual.to_numpy()
        # maxL_inputdata_conceptual = inputsimdata_maxL.l_conceptual.to_numpy()

        # min,max P and L should be the same as in conceptual models
        minP = self.builder.conceptual_model_parameters.get("minP")
        maxP = self.builder.conceptual_model_parameters.get("maxP")
        minL = self.builder.conceptual_model_parameters.get("minL")
        maxL = self.builder.conceptual_model_parameters.get("maxL")

        points_coord_df = pd.DataFrame(points_coord1, columns=["coord1"])
        residual_values = []
        thresholds_values = []
        conceptual_values = []

        # min side
        minL_residual_interpolator = Rbf(
            minL_inputdata_coord1, minL_inputdata_residual, function="linear"
        )
        minL_conceptual_model = conceptual_model(
            points_coord_df, minP=minP, maxP=maxP, minS=minL, maxS=maxL
        )[:, 1]

        minL_minP = np.min(minL_inputdata_coord1)
        minL_minP_val = minL_residual_interpolator(minL_minP)
        minL_maxP = np.max(minL_inputdata_coord1)
        minL_maxP_val = minL_residual_interpolator(minL_maxP)

        residuals = minL_residual_interpolator(points_coord1)
        residuals[points_coord1 > minL_maxP] = minL_maxP_val
        residuals[points_coord1 < minL_minP] = minL_minP_val

        values = minL_conceptual_model - residuals
        values[points_coord1 < minP] = 0
        values[points_coord1 > maxP] = 0

        residual_values.append(residuals)
        thresholds_values.append(values)
        conceptual_values.append(minL_conceptual_model)

        # max side
        maxL_residual_interpolator = Rbf(
            maxL_inputdata_coord1, maxL_inputdata_residual, function="linear"
        )
        maxL_conceptual_model = conceptual_model(
            points_coord_df, minP=minP, maxP=maxP, minS=minL, maxS=maxL
        )[:, 0]

        maxL_minP = np.min(maxL_inputdata_coord1)
        maxL_minP_val = maxL_residual_interpolator(maxL_minP)
        maxL_maxP = np.max(maxL_inputdata_coord1)
        maxL_maxP_val = maxL_residual_interpolator(maxL_maxP)

        residuals = maxL_residual_interpolator(points_coord1)
        residuals[points_coord1 > maxL_maxP] = maxL_maxP_val
        residuals[points_coord1 < maxL_minP] = maxL_minP_val

        values = maxL_conceptual_model - residuals
        values[points_coord1 < minP] = 0
        values[points_coord1 > maxP] = 0

        residual_values.append(residuals)
        thresholds_values.append(values)
        conceptual_values.append(maxL_conceptual_model)

        return thresholds_values, residual_values, conceptual_values

    def interpolate_vertical_thresholds(self, points_coord1, points_coord2):

        function_rbf = "linear"

        conceptual_model = self.builder.vertical_extent_model
        inputsimdata_maxG = self.builder.data_for_vertical_extent_calculation[0]
        inputsimdata_minG = self.builder.data_for_vertical_extent_calculation[1]

        minG_inputdata_coord0 = inputsimdata_minG.coord0.to_numpy()
        minG_inputdata_coord1 = inputsimdata_minG.coord1.to_numpy()
        minG_inputdata_coord2 = inputsimdata_minG.coord2.to_numpy()
        # inputsimdata_minG.coord0.to_numpy()

        # inputsimdata_maxG.coord0.to_numpy()
        maxG_inputdata_coord1 = inputsimdata_maxG.coord1.to_numpy()
        maxG_inputdata_coord2 = inputsimdata_maxG.coord2.to_numpy()
        maxG_inputdata_residual = inputsimdata_maxG.g_residual.to_numpy()
        # inputsimdata_maxG.g_conceptual.to_numpy()

        # min,max P and L should be the same as in conceptual models
        minP = self.builder.conceptual_model_parameters.get("minP")
        maxP = self.builder.conceptual_model_parameters.get("maxP")
        minL = self.builder.conceptual_model_parameters.get("minL")
        maxL = self.builder.conceptual_model_parameters.get("maxL")
        mean_G = self.builder.conceptual_model_parameters.get("mean_growth")
        vertex = self.builder.conceptual_model_parameters.get("vertex")

        points_df = pd.DataFrame()
        points_df["coord1"] = points_coord1
        points_df["coord2"] = points_coord2
        residual_values = []
        threshold_values = []
        conceptual_values = []

        # max growth
        maxG_residual_interpolator = Rbf(
            maxG_inputdata_coord1,
            maxG_inputdata_coord2,
            maxG_inputdata_residual,
            function=function_rbf,
        )
        maxG_conceptual_model = conceptual_model(
            points_df,
            minP=minP,
            maxP=maxP,
            minS=minL,
            maxS=maxL,
            mean_growth=mean_G,
            vertex=vertex,
        )[:, 1]

        # maxG_minP = np.min(maxG_inputdata_coord1)
        # maxG_residual_interpolator(
        #     maxG_minP,
        #     maxG_inputdata_coord2[np.where(maxG_inputdata_coord1 == maxG_minP)][0],
        # )
        # maxG_maxP = np.max(maxG_inputdata_coord1)
        # maxG_residual_interpolator(
        #     maxG_maxP,
        #     maxG_inputdata_coord2[np.where(maxG_inputdata_coord1 == maxG_maxP)][0],
        # )
        # maxG_minL = np.min(maxG_inputdata_coord2)
        # maxG_residual_interpolator(
        #     maxG_inputdata_coord1[np.where(maxG_inputdata_coord2 == maxG_minL)][0],
        #     maxG_minL,
        # )
        # maxG_maxL = np.max(maxG_inputdata_coord2)
        # maxG_residual_interpolator(
        #     maxG_inputdata_coord1[np.where(maxG_inputdata_coord2 == maxG_maxL)][0],
        #     maxG_maxL,
        # )

        residuals = maxG_residual_interpolator(points_coord1, points_coord2)
        thresholds = maxG_conceptual_model - residuals

        residual_values.append(residuals)
        threshold_values.append(thresholds)
        conceptual_values.append(maxG_conceptual_model)

        # intrusion reference contact, conditioning to data
        minG_interpolator = Rbf(
            minG_inputdata_coord1,
            minG_inputdata_coord2,
            minG_inputdata_coord0,
            function=function_rbf,
        )
        thresholds = minG_interpolator(points_coord1, points_coord2)
        minG_conceptual_model = np.zeros(len(points_coord1))

        threshold_values.append(thresholds)
        conceptual_values.append(minG_conceptual_model)

        return threshold_values, residual_values, conceptual_values

    def evaluate_gradient(self, pos):
        ## LG TODO check whether it can be implemented
        raise NotImplementedError("Cannot calculate gradient of Intrusion")

    def evaluate_value(self, pos):
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

        # compute coordinates values for each evaluated point
        intrusion_coord0_pts = self.intrusion_frame[0].evaluate_value(pos)
        intrusion_coord1_pts = self.intrusion_frame[1].evaluate_value(pos)
        intrusion_coord2_pts = self.intrusion_frame[2].evaluate_value(pos)

        self.evaluated_points = [
            pos,
            intrusion_coord0_pts,
            intrusion_coord1_pts,
            intrusion_coord2_pts,
        ]

        thresholds, residuals, conceptual = self.interpolate_lateral_thresholds(
            intrusion_coord1_pts
        )

        if self.intrusion_frame.builder.marginal_faults is not None:
            c2_minside_threshold = thresholds[0]  # np.zeros_like(intrusion_coord2_pts)
            c2_maxside_threshold = thresholds[1]

        else:
            c2_minside_threshold = thresholds[0]
            c2_maxside_threshold = thresholds[1]

        thresholds, residuals, conceptual = self.interpolate_vertical_thresholds(
            intrusion_coord1_pts, intrusion_coord2_pts
        )
        c0_minside_threshold = thresholds[1]
        c0_maxside_threshold = thresholds[0]

        if len(self.assisting_faults) > 0:
            fault = self.assisting_faults.get("structure")
            weight = self.assisting_faults.get("asymmetry_weight", 1)
            evaluation_points_in_fault = fault[0].evaluate_value(pos)
            c0_maxside_threshold[evaluation_points_in_fault >= 0] = (
                c0_maxside_threshold[evaluation_points_in_fault >= 0] * weight
            )

        mid_point = c0_minside_threshold + ((c0_maxside_threshold - c0_minside_threshold) / 2)

        a = intrusion_coord2_pts >= c2_maxside_threshold
        b = intrusion_coord2_pts <= c2_minside_threshold
        c = (
            (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
            * (intrusion_coord0_pts <= c0_minside_threshold)
        )
        d = (
            (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
            * (intrusion_coord0_pts >= c0_maxside_threshold)
        )
        e = (
            (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
            * (mid_point >= intrusion_coord0_pts)
            * (intrusion_coord0_pts > c0_minside_threshold)
        )
        f = (
            (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
            * (mid_point < intrusion_coord0_pts)
            * (intrusion_coord0_pts < c0_maxside_threshold)
        )

        mod_Smin_thresholds = intrusion_coord2_pts - c2_minside_threshold
        mod_Smax_thresholds = intrusion_coord2_pts - c2_maxside_threshold
        mod_Gmin_thresholds = intrusion_coord0_pts - c0_minside_threshold
        mod_Gmax_thresholds = intrusion_coord0_pts - c0_maxside_threshold

        intrusion_sf = (
            a * mod_Smax_thresholds
            + b * abs(mod_Smin_thresholds)
            + c * abs(mod_Gmin_thresholds)
            + d * mod_Gmax_thresholds
            - e * mod_Gmin_thresholds
            + f * mod_Gmax_thresholds
        ) * (
            -1
        )  # multiply by (-1) so intrusions can be used as unconformities

        return intrusion_sf

    def evaluate_value_test(self, points):
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

        # compute coordinates values for each evaluated point
        intrusion_coord0_pts = self.intrusion_frame[0].evaluate_value(points)
        intrusion_coord1_pts = self.intrusion_frame[1].evaluate_value(points)
        intrusion_coord2_pts = self.intrusion_frame[2].evaluate_value(points)

        self.evaluated_points = [
            points,
            intrusion_coord0_pts,
            intrusion_coord1_pts,
            intrusion_coord2_pts,
        ]

        thresholds, residuals, conceptual = self.interpolate_lateral_thresholds(
            intrusion_coord1_pts
        )

        if self.intrusion_frame.builder.marginal_faults is not None:
            c2_minside_threshold = np.zeros_like(intrusion_coord2_pts)
            c2_maxside_threshold = thresholds[1]

        else:
            c2_minside_threshold = thresholds[0]
            c2_maxside_threshold = thresholds[1]

        thresholds, residuals, conceptual = self.interpolate_vertical_thresholds(
            intrusion_coord1_pts, intrusion_coord2_pts
        )
        c0_minside_threshold = thresholds[1]
        c0_maxside_threshold = thresholds[0]

        mid_point = c0_minside_threshold + ((c0_maxside_threshold - c0_minside_threshold) / 2)

        mod_intrusion_coord0_pts = intrusion_coord0_pts - mid_point
        mod_c0_minside_threshold = c0_minside_threshold - mid_point
        mod_c0_maxside_threshold = c0_maxside_threshold + mid_point

        a = (
            (mod_intrusion_coord0_pts >= mid_point)
            * (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
        )
        b = (
            (mod_intrusion_coord0_pts <= mid_point)
            * (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
        )
        c = (
            (mod_intrusion_coord0_pts <= mid_point)
            * (mod_intrusion_coord0_pts >= mod_c0_minside_threshold)
            * (c2_minside_threshold < intrusion_coord2_pts)
            * (intrusion_coord2_pts < c2_maxside_threshold)
        )

        intrusion_sf = mod_intrusion_coord0_pts
        intrusion_sf[a] = mod_intrusion_coord0_pts[a] - mod_c0_maxside_threshold[a]
        intrusion_sf[b] = abs(mod_c0_minside_threshold[b] + mod_intrusion_coord0_pts[b])
        intrusion_sf[c] = mod_intrusion_coord0_pts[c] - mod_c0_minside_threshold[c]

        return intrusion_sf

    def get_data(self, value_map: Optional[dict] = None):
        pass

    def copy(self):
        pass
