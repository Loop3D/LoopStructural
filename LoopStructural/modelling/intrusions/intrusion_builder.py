import numpy as np
import pandas as pd

from ...utils import getLogger
from .intrusion_feature import IntrusionFeature


from ..features.builders import BaseBuilder
from ...utils import rng
from .geometric_scaling_functions import *

logger = getLogger(__name__)


class IntrusionBuilder(BaseBuilder):
    def __init__(
        self,
        frame,
        model,
        vertical_extent_model=None,
        lateral_extent_model=None,
        name="intrusion builder",
        **kwargs,
    ):
        """
        Constructor for an IntrusionBuilder
        Sets up interpolators of vertical and lateral contacts

        Parameters
        ----------
        frame : StructuralFrame
            the structural frame that the intrusion will be built using
        model : GeologicalModel, optional
            the geological model that the intrusion will belong to, by default None
        name : str, optional
            name of the intrusion, by default "intrusion builder"
        interpolator : GeologicalInterpolator, optional
            a geological interpolator, by default None
        """

        BaseBuilder.__init__(self, model, name=name)

        self.intrusion_frame = frame
        self._up_to_date = False
        self._feature = IntrusionFeature(
            frame=frame,
            builder=self,
            name=self.name,
        )

        self._build_arguments = {}
        self.data = None
        self.data_prepared = False
        self.lateral_contact_data = None
        self.vertical_contact_data = None
        self.lateral_extent_model = lateral_extent_model
        self.vertical_extent_model = vertical_extent_model
        self.width_data = [True, True]
        self.thickness_data = False
        self.constrain_sides_with_rooffloor_data = False

        self.data_for_lateral_extent_calculation = None
        self.data_for_vertical_extent_calculation = None
        self.evaluation_grid = None
        self.conceptual_model_parameters = {}

        self.marginal_faults = self._feature.intrusion_frame.builder.marginal_faults

    def create_grid_for_evaluation(self, spacing=None):
        """
        Create the grid points in which to simulate vertical and lateral

        Parameters
        ----------
        spacing : np.array
            list/array with spacing value for X,Y,Z

        Returns
        -------
        """

        if spacing is None:
            spacing = self.model.nsteps

        grid_points = self.model.regular_grid(nsteps=spacing, shuffle=False)

        grid_points_coord0 = self.intrusion_frame[0].evaluate_value(grid_points)

        grid_points_coord1 = self.intrusion_frame[1].evaluate_value(grid_points)

        grid_points_coord2 = self.intrusion_frame[2].evaluate_value(grid_points)

        self.evaluation_grid = [
            grid_points,
            grid_points_coord0,
            grid_points_coord1,
            grid_points_coord2,
            spacing,
        ]

    def set_data_for_extent_calculation(self, intrusion_data: pd.DataFrame):
        """Set data for lateral extent (distances in c2 axis)  and vertical extent (distances in c0 axis) simulation.
        creates a copy of the data

        Parameters
        ----------
        intrusion_data : DataFrame

        Returns
        ----------

        """

        self.data = intrusion_data.copy()

    def create_geometry_using_geometric_scaling(
        self, geometric_scaling_parameters, reference_contact_data
    ):

        geometric_scaling_parameters.get("intrusion_type", None)
        intrusion_length = geometric_scaling_parameters.get("intrusion_length", None)
        geometric_scaling_parameters.get("inflation_vector", np.array([[0, 0, 1]]))
        thickness = geometric_scaling_parameters.get("thickness", None)

        if (
            self.intrusion_frame.builder.intrusion_network_contact == "floor"
            or self.intrusion_frame.builder.intrusion_network_contact == "base"
        ):
            geometric_scaling_parameters.get("inflation_vector", np.array([[0, 0, 1]]))
        else:
            geometric_scaling_parameters.get("inflation_vector", np.array([[0, 0, -1]]))

        if intrusion_length is None and thickness is None:
            raise ValueError(
                "No {} data. Add intrusion_type and intrusion_length (or thickness) to geometric_scaling_parameters dictionary".format(
                    self.intrusion_frame.builder.intrusion_other_contact
                )
            )

        else:  # -- create synthetic data to constrain interpolation using geometric scaling
            estimated_thickness = thickness
            if estimated_thickness is None:
                raise Exception('Not implemented')
                # estimated_thickness = thickness_from_geometric_scaling(
                #     intrusion_length, intrusion_type
                # )

            print(
                "Building tabular intrusion using geometric scaling parameters: estimated thicknes = {} meters".format(
                    round(estimated_thickness)
                )
            )
            raise Exception('Not implemented')
            # (
            #     other_contact_data_temp,
            #     other_contact_data_xyz_temp,
            # ) = contact_pts_using_geometric_scaling(
            #     estimated_thickness, reference_contact_data, inflation_vector
            # )

            # return other_contact_data_temp

    def prepare_data(self, geometric_scaling_parameters):
        """Prepare the data to compute distance thresholds along the frame coordinates.
        These distance thresholds represent the contact of the intrusion.

        1. Select lateral data and separate it between sides (i.e., c2>0 or c2<0)
        2. Separate data between roof and floor contact data.
            If there are only data for one of the contact (only roof or floor contact),
            it can use geometric scaling to create synthetic data of the opposite contact.
            This assumes constant thickness.

        """
        if self.data is None or self.data.shape[0] == 0:
            raise ValueError("Cannot create intrusion with no data")

        intrusion_data = self.data.copy()

        data_xyz = intrusion_data.loc[:, ["X", "Y", "Z"]].to_numpy()
        intrusion_data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(data_xyz)
        intrusion_data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(data_xyz)
        intrusion_data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(data_xyz)

        # -- separate data between both sides of the intrusion, using intrusion axis (i.e., coord2 = 0)

        data_minside = intrusion_data[
            (intrusion_data["intrusion_side"]) & (intrusion_data["coord2"] <= 0)
        ].copy()
        data_minside.reset_index(inplace=True, drop=True)

        data_maxside = intrusion_data[
            (intrusion_data["intrusion_side"]) & (intrusion_data["coord2"] > 0)
        ].copy()
        data_maxside.reset_index(inplace=True, drop=True)

        if data_minside.shape[0] < 3:  # minimum three input data for interpolation
            self.width_data[0] = False

        else:
            self.width_data[0] = True

        if data_maxside.shape[0] < 3:  # minimum three input data for interpolation
            self.width_data[1] = False

        else:
            self.width_data[1] = True

        data_sides = pd.concat([data_minside, data_maxside])
        data_sides.reset_index(inplace=True, drop=True)

        self.lateral_contact_data = [data_sides, data_minside, data_maxside]

        # -- separate data between roof and floor data

        intrusion_network_data_xyz = self.intrusion_frame.builder.intrusion_network_data.loc[
            :, ["X", "Y", "Z"]
        ].to_numpy()
        intrusion_network_data = self.intrusion_frame.builder.intrusion_network_data.loc[
            :, ["X", "Y", "Z"]
        ].copy()
        intrusion_network_data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(
            intrusion_network_data_xyz
        )
        intrusion_network_data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(
            intrusion_network_data_xyz
        )
        intrusion_network_data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(
            intrusion_network_data_xyz
        )
        intrusion_network_data.reset_index(inplace=True)

        # -- if no data points for roof or floor, use geometric scaling to create points for SGS
        if self.intrusion_frame.builder.other_contact_data.shape[0] == 0:
            if len(geometric_scaling_parameters) == 0:
                # self.create_geometry_from_conceptual_model()
                self.thickness_data = False  # set to False, so other_contact is constrained only with conceptual model.
                other_contact_data_xyz = intrusion_network_data_xyz
                other_contact_data = intrusion_network_data

            else:
                other_contact_data_temp1 = self.intrusion_frame.builder.other_contact_data
                other_contact_data_temp2 = self.create_geometry_using_geometric_scaling(
                    geometric_scaling_parameters, intrusion_network_data
                )

                other_contact_data = pd.concat([other_contact_data_temp1, other_contact_data_temp2])

                other_contact_data_xyz = other_contact_data.loc[:, ["X", "Y", "Z"]].to_numpy()

        else:
            self.thickness_data = True
            other_contact_data_xyz = self.intrusion_frame.builder.other_contact_data.loc[
                :, ["X", "Y", "Z"]
            ].to_numpy()
            other_contact_data = self.intrusion_frame.builder.other_contact_data.loc[
                :, ["X", "Y", "Z"]
            ].copy()

        other_contact_data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(
            other_contact_data_xyz
        )
        other_contact_data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(
            other_contact_data_xyz
        )
        other_contact_data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(
            other_contact_data_xyz
        )
        other_contact_data.reset_index(inplace=True)

        self.vertical_contact_data = [intrusion_network_data, other_contact_data]
        self.data_prepared = True

    def set_conceptual_models_parameters(self):
        """Creates a dictionary of parameters used for conceptual models.
        These parameters includes the basic parameters for the functions
        representing the conceptual models of the intrusion geometry.

        """
        if not callable(self.lateral_extent_model) or not callable(self.vertical_extent_model):
            raise ValueError("lateral_extent_model and vertical_extent_model must be functions")

        grid_points_coord1 = self.evaluation_grid[2]

        modelcover, minP, maxP, minL, maxL = self.lateral_extent_model()
        mean_c0 = self.vertical_extent_model()

        if minL is None:
            minL = min(
                self.vertical_contact_data[0]["coord2"].min(),
                self.vertical_contact_data[1]["coord2"].min(),
                self.lateral_contact_data[0]["coord2"].min(),
            )

        if maxL is None:
            maxL = max(
                self.vertical_contact_data[0]["coord2"].max(),
                self.vertical_contact_data[1]["coord2"].max(),
                self.lateral_contact_data[0]["coord2"].max(),
            )

        if minL < 0 and maxL < 0:
            maxL = minL * -1

        if minL > 0 and maxL > 0:
            minL = maxL * -1

        if modelcover is True:
            minP = np.nanmin(grid_points_coord1)
            maxP = np.nanmax(grid_points_coord1)
        else:
            if minP is None:
                minP = min(
                    self.vertical_contact_data[0]["coord1"].min(),
                    self.vertical_contact_data[1]["coord1"].min(),
                    self.lateral_contact_data[0]["coord1"].min(),
                )
            if maxP is None:
                maxP = max(
                    self.vertical_contact_data[0]["coord1"].max(),
                    self.vertical_contact_data[1]["coord1"].max(),
                    self.lateral_contact_data[0]["coord1"].max(),
                )

        # extra parameters for growth
        if mean_c0 is None:
            mean_growth = self.vertical_contact_data[1].loc[:, "coord0"].mean()
        else:
            mean_growth = mean_c0

        maxG = self.vertical_contact_data[1]["coord0"].max()
        coord_PL_for_maxG = (
            self.vertical_contact_data[1][
                self.vertical_contact_data[1].coord0 == self.vertical_contact_data[1].coord0.max()
            ]
            .loc[:, ["coord1", "coord2"]]
            .to_numpy()
        )

        self.conceptual_model_parameters["minP"] = minP
        self.conceptual_model_parameters["maxP"] = maxP
        self.conceptual_model_parameters["minL"] = minL
        self.conceptual_model_parameters["maxL"] = maxL
        self.conceptual_model_parameters["model_cover"] = modelcover
        self.conceptual_model_parameters["mean_growth"] = mean_growth
        self.conceptual_model_parameters["vertex"] = [
            coord_PL_for_maxG[0][0],
            coord_PL_for_maxG[0][1],
            maxG,
        ]

    def set_data_for_lateral_thresholds(self):
        """
        Sets the data for the interpolation of lateral thresholds.
        Computes conceptual and residual values

        Parameters
        ----------

        Returns
        -------
        """

        # -- generate data frame containing input data for simulation

        self.set_conceptual_models_parameters()

        minP = self.conceptual_model_parameters.get("minP")
        maxP = self.conceptual_model_parameters.get("maxP")
        minL = self.conceptual_model_parameters.get("minL")
        maxL = self.conceptual_model_parameters.get("maxL")

        if self.width_data[0] is False:  # i.e., no lateral data for side L<0
            print(
                "Not enought lateral data to constrain side L<0. Conceptual model will be used to constrain lateral extent"
            )

            random_p = pd.DataFrame(rng.uniform(minP, maxP, 10), columns=["coord1"])
            conceptual_l = self.lateral_extent_model(
                lateral_contact_data=random_p,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )
            data_for_min_L = pd.DataFrame(
                np.vstack([conceptual_l[:, 1], random_p.loc[:, "coord1"].to_numpy()]).T,
                columns=["l_conceptual", "coord1"],
            )
            data_for_min_L.loc[:, "l_residual"] = 0

            if len(self.lateral_contact_data[1]) > 0:
                data_minL = self.lateral_contact_data[1]
                data_conceptual_minL = self.lateral_extent_model(
                    lateral_contact_data=data_minL,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                data_residual_minL = (
                    data_conceptual_minL[:, 1] - data_minL.loc[:, "coord2"]
                ).to_numpy()
                data_for_min_L_ = data_minL.loc[
                    :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
                ].copy()
                data_for_min_L_.loc[:, "l_residual"] = data_residual_minL
                data_for_min_L_.loc[:, "l_conceptual"] = data_conceptual_minL[:, 1]
                data_for_min_L_.reset_index(inplace=True)
                data_for_min_L = pd.concat([data_for_min_L, data_for_min_L_])

            data_for_min_L["l_residual"] = data_for_min_L["l_residual"].astype(float)
            data_for_min_L["coord1"] = data_for_min_L["coord1"].astype(float)

        else:
            # -- Side of intrusion with coord2<0 (l<0)
            data_minL = self.lateral_contact_data[1]
            data_conceptual_minL = self.lateral_extent_model(
                lateral_contact_data=data_minL,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )
            data_residual_minL = (
                data_conceptual_minL[:, 1] - data_minL.loc[:, "coord2"]
            ).to_numpy()
            data_for_min_L = data_minL.loc[:, ["X", "Y", "Z", "coord0", "coord1", "coord2"]].copy()
            data_for_min_L.loc[:, "l_residual"] = data_residual_minL
            data_for_min_L.loc[:, "l_conceptual"] = data_conceptual_minL[:, 1]
            data_for_min_L.reset_index(inplace=True)
            # data_for_min_L.loc[:, "ref_coord"] = 0

        if not self.width_data[1]:  # i.e., no lateral data for side L>0
            print(
                "Not enought lateral data to constrain side L>0. Conceptual model will be used to constrain lateral extent"
            )

            random_p = pd.DataFrame(rng.uniform(minP, maxP, 10), columns=["coord1"])
            conceptual_l = self.lateral_extent_model(
                lateral_contact_data=random_p,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )
            data_for_max_L = pd.DataFrame(
                np.vstack([conceptual_l[:, 0], random_p.loc[:, "coord1"].to_numpy()]).T,
                columns=["l_conceptual", "coord1"],
            )
            data_for_max_L.loc[:, "l_residual"] = 0

            if len(self.lateral_contact_data[2]) > 0:
                data_maxL = self.lateral_contact_data[2]
                data_conceptual_maxL = self.lateral_extent_model(
                    lateral_contact_data=data_maxL,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                data_residual_maxL = (
                    data_conceptual_maxL[:, 0] - data_maxL.loc[:, "coord2"]
                ).to_numpy()
                data_for_max_L_ = data_maxL.loc[
                    :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
                ].copy()
                data_for_max_L_.loc[:, "l_residual"] = data_residual_maxL
                data_for_max_L_.loc[:, "l_conceptual"] = data_conceptual_maxL[:, 0]
                data_for_max_L_.reset_index(inplace=True)
                data_for_max_L = pd.concat([data_for_max_L, data_for_max_L_])

            data_for_max_L["l_residual"] = data_for_max_L["l_residual"].astype(float)
            data_for_max_L["coord1"] = data_for_max_L["coord1"].astype(float)

        else:
            data_maxL = self.lateral_contact_data[2]
            data_conceptual_maxL = self.lateral_extent_model(
                lateral_contact_data=data_maxL,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )
            data_residual_maxL = (
                data_conceptual_maxL[:, 0] - data_maxL.loc[:, "coord2"]
            ).to_numpy()
            data_for_max_L = data_maxL.loc[:, ["X", "Y", "Z", "coord0", "coord1", "coord2"]].copy()
            data_for_max_L.loc[:, "l_residual"] = data_residual_maxL
            data_for_max_L.loc[:, "l_conceptual"] = data_conceptual_maxL[:, 0]
            data_for_max_L.reset_index(inplace=True)
            # data_for_max_L.loc[:, "ref_coord"] = 0

        # check if roof or floor data outside of conceptual model.
        # if so, add as constraints to conceptual model.

        vertical_data = pd.concat([self.vertical_contact_data[0], self.vertical_contact_data[1]])
        vertical_data.loc[:, ["conceptual_maxside", "conceptual_minside"]] = (
            self.lateral_extent_model(
                lateral_contact_data=vertical_data,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )
        )

        data_minL_temp = vertical_data[vertical_data["coord2"] < 0].copy()
        data_for_min_L_ = (
            data_minL_temp[data_minL_temp["coord2"] < data_minL_temp["conceptual_minside"]]
            .loc[:, ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_minside"]]
            .copy()
        )
        data_for_min_L_.loc[:, "l_residual"] = (
            data_for_min_L_.loc[:, "conceptual_minside"] - data_for_min_L_.loc[:, "coord2"]
        )
        data_for_min_L_.rename(columns={"conceptual_minside": "l_conceptual"}, inplace=True)
        data_for_min_L_.reset_index(inplace=True)
        data_for_min_L_.drop_duplicates(
            subset=[
                "X",
                "Y",
                "Z",
                "coord0",
                "coord1",
                "coord2",
                "l_conceptual",
                "l_residual",
            ],
            inplace=True,
        )

        if len(data_for_min_L_) > 0 and self.constrain_sides_with_rooffloor_data:
            print("adding data from roof/floor to constrain L<0")
            data_for_min_L = pd.concat([data_for_min_L, data_for_min_L_])

        data_maxL_temp = vertical_data[vertical_data["coord2"] >= 0].copy()
        data_for_max_L_ = (
            data_maxL_temp[data_maxL_temp["coord2"] > data_maxL_temp["conceptual_maxside"]]
            .loc[:, ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_maxside"]]
            .copy()
        )
        data_for_max_L_.loc[:, "l_residual"] = (
            data_for_max_L_.loc[:, "conceptual_maxside"] - data_for_max_L_.loc[:, "coord2"]
        )
        data_for_max_L_.rename(columns={"conceptual_maxside": "l_conceptual"}, inplace=True)
        data_for_max_L_.reset_index(inplace=True)
        data_for_max_L_.drop_duplicates(
            subset=[
                "X",
                "Y",
                "Z",
                "coord0",
                "coord1",
                "coord2",
                "l_conceptual",
                "l_residual",
            ],
            inplace=True,
        )

        if len(data_for_max_L_) > 0 and self.constrain_sides_with_rooffloor_data:
            print("adding data from roof/floor to constrain L>0")
            data_for_max_L = pd.concat([data_for_max_L, data_for_max_L_])

        data_for_min_L["l_residual"] = data_for_min_L["l_residual"].astype(float)
        data_for_min_L["coord1"] = data_for_min_L["coord1"].astype(float)
        data_for_max_L["l_residual"] = data_for_max_L["l_residual"].astype(float)
        data_for_max_L["coord1"] = data_for_max_L["coord1"].astype(float)

        self.data_for_lateral_extent_calculation = [data_for_min_L, data_for_max_L]

    def set_data_for_vertical_thresholds(self):
        """
        Sets the data for the interpolation of the vertical thresholds.
        Computes conceptual and residual values

        Parameters
        ----------

        Returns
        -------
        """

        # -- Generate data frame containing input data for simulation
        inet_data = self.vertical_contact_data[0]
        other_contact_data = self.vertical_contact_data[1]

        # # --- parameters for conceptual model
        minP = self.conceptual_model_parameters.get("minP")
        maxP = self.conceptual_model_parameters.get("maxP")
        minL = self.conceptual_model_parameters.get("minL")
        maxL = self.conceptual_model_parameters.get("maxL")
        meanG = self.conceptual_model_parameters.get("mean_growth")
        vertex = self.conceptual_model_parameters.get("vertex")

        # --- growth simulation input data (max G, simulation of contact opposite to intrusion network)

        data_conceptual_G = self.vertical_extent_model(
            other_contact_data,
            mean_growth=meanG,
            minP=minP,
            maxP=maxP,
            minS=minL,
            maxS=maxL,
            vertex=vertex,
        )
        data_residual_G = (data_conceptual_G[:, 1] - other_contact_data.loc[:, "coord0"]).to_numpy()
        inputsimdata_maxG = other_contact_data.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()
        inputsimdata_maxG.loc[:, "g_residual"] = data_residual_G
        inputsimdata_maxG.loc[:, "g_conceptual"] = data_conceptual_G[:, 1]

        if not self.thickness_data:
            inputsimdata_maxG.loc[:, "g_residual"] = 0

            p = np.linspace(minP, maxP, 10)
            l = np.linspace(minL, maxL, 10)

            pp, ll = np.meshgrid(p, l)
            pl = np.array([pp.flatten(), ll.flatten()]).T
            rng.shuffle(pl)

            inputsimdata_maxG_ = pd.DataFrame(pl[:30, :], columns=["coord1", "coord2"])
            data_conceptual_G_ = self.vertical_extent_model(
                inputsimdata_maxG_,
                mean_growth=meanG,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
                vertex=vertex,
            )

            inputsimdata_maxG_.loc[:, "g_conceptual"] = data_conceptual_G_[:, 1]
            inputsimdata_maxG_.loc[:, "g_residual"] = 0

            inputsimdata_maxG_complete = pd.concat([inputsimdata_maxG, inputsimdata_maxG_])

        else:
            inputsimdata_maxG_complete = inputsimdata_maxG

        # --- growth simulation input data for intrusion network conditioning
        inputsimdata_inetG = inet_data.loc[:, ["X", "Y", "Z", "coord0", "coord1", "coord2"]].copy()

        self.data_for_vertical_extent_calculation = [
            inputsimdata_maxG_complete,
            inputsimdata_inetG,
        ]

    def build(
        self,
        # parameters_for_extent_sgs={},
        geometric_scaling_parameters={},
        **kwargs,
    ):
        """Main building function for intrusion.
        Set up interpolators for extent calculation
        If SGS --> Calculates variogram and simulate thresholds along frame axes

        Parameters
        ----------
        vertical_extent_sgs_parameters : dict, optional
            parameters for the vertical sequential gaussian simulation, by default {}
        lateral_extent_sgs_parameters : dict, optional
            parameters for the vertical sequential gaussian simulation, by default {}
        """
        self.prepare_data(geometric_scaling_parameters)
        self.create_grid_for_evaluation()

        self.set_data_for_lateral_thresholds()
        self.set_data_for_vertical_thresholds()
