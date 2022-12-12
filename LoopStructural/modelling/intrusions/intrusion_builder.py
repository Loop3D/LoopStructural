import numpy as np
import pandas as pd

from ...utils import getLogger
from .intrusion_feature import IntrusionFeature
from ...interpolators import StructuredGrid2D

from scipy.interpolate import Rbf

from ...modelling.intrusions.intrusion_support_functions import (
    grid_from_array,
)

from .geometric_scaling_functions import *

logger = getLogger(__name__)

# import GSLIB library
try:
    import geostatspy.GSLIB as GSLIB  # GSLIB utilities, viz and wrapped functions

except ImportError:
    logger.warning(
        "Cannot use Intrusions: GeostatPy not installed \n" "pip install geostatspy"
    )
    raise ImportError("GSLIB")
try:
    import geostatspy.geostats as geostats  # GSLIB converted to Python
except ImportError:
    logger.warning(
        "Cannot use Intrusions: GeostatPy not installed \n" "pip install geostatspy"
    )
    raise ImportError("geostats")


class IntrusionBuilder:
    def __init__(self, frame, model=None, name="intrusion"):
        self.name = name
        self.intrusion_frame = frame
        self._up_to_date = False
        self.faults = []
        self._feature = IntrusionFeature(
            frame=frame, builder=self, faults=self.faults, name=self.name
        )

        # contact data:
        self.lateral_contact_data = None
        self.vertical_contact_data = None

        # sequential gaussian simulation parameters:
        self.lateral_sgs_parameters = None
        self.vertical_sgs_parameters = None

        self.lateral_sgs_variogram = None
        self.vertical_sgs_variogram = None

        self.lateral_sgs_input_data = None
        self.vertical_sgs_input_data = None

        # self.lateral_extent_sgs_parameters = None
        # self.vertical_extent_sgs_parameters = None

        self.simulation_grid = None
        self.data = None
        self.data_prepared = False
        self.model = model
        self._build_arguments = {}
        self.width_data = [True, True]

        # self.simulationGSLIB_s_outcome = None
        self.growth_simulated_thresholds_grid = None
        self.conceptual_model_parameters = {}

        self.lateral_conceptual_model_parameters = {}

    @property
    def feature(self):
        return self._feature

    @property
    def build_arguments(self):
        return self._build_arguments

    @build_arguments.setter
    def build_arguments(self, arguments):
        """Set the build arguments and flag that
        up to date is False

        Parameters
        ----------
        arguments : dictionary
            dictionary containing keys for variogram arguments
        """
        if type(arguments) == dict:
            self._up_to_date = False
            self._build_arguments = arguments
        else:
            logger.error(
                f"Cannot update build arguments with {type(arguments)}, must be a dictionary"
            )

    def add_fault(self, fault):
        """
        Add a fault to the intrusion feature builder

        Parameters
        ----------
        fault : FaultSegment
            A faultsegment to add to the geological feature

        Returns
        -------

        """
        self._up_to_date = False
        self.faults.append(fault)

    def create_grid_for_simulation(self, spacing=None):
        """
        Create the grid points in which to simulate vertical and lateral

        Parameters
        ----------
        spacing = list/array with spacing value for X,Y,Z

        Returns
        -------
        """

        if spacing == None:
            spacing = self.model.nsteps

        grid_points = self.model.regular_grid(spacing, shuffle=False)

        grid_points_coord0 = self.intrusion_frame[0].evaluate_value(grid_points)

        grid_points_coord1 = self.intrusion_frame[1].evaluate_value(grid_points)

        grid_points_coord2 = self.intrusion_frame[2].evaluate_value(grid_points)

        self.simulation_grid = [
            grid_points,
            grid_points_coord0,
            grid_points_coord1,
            grid_points_coord2,
            spacing,
        ]

    def set_data_for_extent_simulation(self, intrusion_data):
        """Set data for lateral extent (distances in L axis)  and vertical extent (distances in G axis) simulation.

        Parameters
        ----------
        intrusion_data : DataFrame

        Returns
        ----------

        """

        self.data = intrusion_data.copy()

    def prepare_data(self, geometric_scaling_parameters):
        """ """
        if self.data is None or self.data.shape[0] == 0:
            raise ValueError("Cannot create intrusion with no data")
        data_xyz = self.data.loc[:, ["X", "Y", "Z"]].to_numpy()
        self.data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(data_xyz)
        self.data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(data_xyz)
        self.data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(data_xyz)

        # -- separate data between both sides of the intrusion, using intrusion axis (i.e. coord2 = 0)

        data_minside = self.data[
            (self.data["intrusion_side"] == True) & (self.data["coord2"] <= 0)
        ].copy()
        data_minside.reset_index(inplace=True, drop=True)

        data_maxside = self.data[
            (self.data["intrusion_side"] == True) & (self.data["coord2"] > 0)
        ].copy()
        data_maxside.reset_index(inplace=True, drop=True)

        if data_minside.shape[0] < 3:  # minimum three inut data for SGS simulation
            self.width_data[0] = False

        else:
            self.width_data[0] = True

        if data_maxside.shape[0] < 3:  # minimum three inut data for SGS simulation
            self.width_data[1] = False

        else:
            self.width_data[1] = True

        data_sides = pd.concat([data_minside, data_maxside])
        data_sides.reset_index(inplace=True, drop=True)

        self.lateral_contact_data = [data_sides, data_minside, data_maxside]

        # -- separate data between roof and floor data

        intrusion_network_data_xyz = (
            self.intrusion_frame.builder.intrusion_network_data.loc[
                :, ["X", "Y", "Z"]
            ].to_numpy()
        )
        intrusion_network_data = (
            self.intrusion_frame.builder.intrusion_network_data.loc[
                :, ["X", "Y", "Z"]
            ].copy()
        )
        intrusion_network_data.loc[:, "coord0"] = self.intrusion_frame[
            0
        ].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:, "coord1"] = self.intrusion_frame[
            1
        ].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:, "coord2"] = self.intrusion_frame[
            2
        ].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.reset_index(inplace=True)

        # -- if no data points for roof or floor, use geometric scaling to create points for SGS
        if self.intrusion_frame.builder.other_contact_data.shape[0] < 3:
            other_contact_data_temp1 = self.intrusion_frame.builder.other_contact_data

            intrusion_type = geometric_scaling_parameters.get("intrusion_type", None)
            intrusion_length = geometric_scaling_parameters.get(
                "intrusion_length", None
            )
            inflation_vector = geometric_scaling_parameters.get(
                "inflation_vector", np.array([[0, 0, 1]])
            )
            thickness = geometric_scaling_parameters.get("thickness", None)

            if (
                self.intrusion_frame.builder.intrusion_network_contact == "floor"
                or self.intrusion_frame.builder.intrusion_network_contact == "base"
            ):
                inflation_vector = geometric_scaling_parameters.get(
                    "inflation_vector", np.array([[0, 0, 1]])
                )
            else:
                inflation_vector = geometric_scaling_parameters.get(
                    "inflation_vector", np.array([[0, 0, -1]])
                )

            if intrusion_length == None and thickness == None:
                raise ValueError(
                    "No {} data. Add intrusion_type and intrusion_length (or thickness) to geometric_scaling_parameters dictionary".format(
                        self.intrusion_frame.builder.intrusion_other_contact
                    )
                )

            else:  # -- create data using geometric scaling
                estimated_thickness = thickness
                if estimated_thickness == None:
                    estimated_thickness = thickness_from_geometric_scaling(
                        intrusion_length, intrusion_type
                    )

                print(
                    "Building intrusion using geometric scaling parameters: estimated thicknes = {} meters".format(
                        round(estimated_thickness)
                    )
                )
                (
                    other_contact_data_temp2,
                    other_contact_data_xyz_temp,
                ) = contact_pts_using_geometric_scaling(
                    estimated_thickness, intrusion_network_data, inflation_vector
                )

                other_contact_data = pd.concat(
                    [other_contact_data_temp1, other_contact_data_temp2]
                )
                other_contact_data_xyz = other_contact_data.loc[
                    :, ["X", "Y", "Z"]
                ].to_numpy()

        else:
            other_contact_data_xyz = (
                self.intrusion_frame.builder.other_contact_data.loc[
                    :, ["X", "Y", "Z"]
                ].to_numpy()
            )
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

    def set_l_sgs_GSLIBparameters(self, lateral_simulation_parameters):
        """
        Simulation parameters for lateral extent simulation

        Parameters
        ----------

        lateral_simulation_parameters: python dictionary with the following parameters -->

        tmin, tmax : all values, regardless of which variable, strictly less than tmin and greater than or equal to tmax are ignored.
        itrans : 0 - no transformation requiered, data already nscored / 1 - for transformation
        ktype : type of interpolation, 0 for simple kriging - 1 for ordinary kriging

        nx, ny : Numbers of blocks. Grid node indices ix, iy increase from 1 to nx, ny respectively (in the positive x,y direction).
        xmn, ymn :  Centre of the first block (xmn,ymn)
        xsiz, ysiz : Size of blocks
        zmin, zmax : The minimum and maximum allowable data values. These are used in the back-transformation procedure.
        nxdis, nydis : Number of discretization points for a block. If both nxdis and nydis are set to 1, then point kriging is performed.
        ndmin, ndmax : The minimum and maximum number of original data that should be used to simulate a grid node. If there are fewer than ndmin data points, the node is not simulated.
        radius : The maximum isotropic search radius.

        Returns
        ------------
        """
        self._up_to_date = False
        tmin = lateral_simulation_parameters.get("tmin", -9999)
        tmax = lateral_simulation_parameters.get("tmax", 9999)
        itrans = lateral_simulation_parameters.get("itrans", 1)
        ktype = lateral_simulation_parameters.get("ktype", 0)
        nx = lateral_simulation_parameters.get("nx", 200)
        ny = lateral_simulation_parameters.get("ny", 3)
        xmn = lateral_simulation_parameters.get("xmn", None)
        ymn = lateral_simulation_parameters.get("ymn", -0.0001)
        xsiz = lateral_simulation_parameters.get("xsiz", None)
        ysiz = lateral_simulation_parameters.get("ysiz", 0.0001)
        zmin = lateral_simulation_parameters.get("zmin", None)
        zmax = lateral_simulation_parameters.get("zmax", None)
        nxdis = lateral_simulation_parameters.get("nxdis", 1)
        nydis = lateral_simulation_parameters.get("nydis", 1)
        ndmin = lateral_simulation_parameters.get("ndmin", 0)
        ndmax = lateral_simulation_parameters.get("ndmax", 3)
        radius = lateral_simulation_parameters.get("radius", 500)
        nugget = lateral_simulation_parameters.get("nugget", 0)
        nst = lateral_simulation_parameters.get("nst", 1)
        it1 = lateral_simulation_parameters.get("it1", 1)
        cc1 = lateral_simulation_parameters.get("cc1", 1)
        azi1 = lateral_simulation_parameters.get("azi1", 90)
        hmaj1 = lateral_simulation_parameters.get("hmaj1", 999999)
        hmin1 = lateral_simulation_parameters.get("hmin1", 999999)

        l_parameters = {
            "tmin": tmin,
            "tmax": tmax,
            "itrans": itrans,
            "ktype": ktype,
            "nx": nx,
            "ny": ny,
            "xmn": xmn,
            "ymn": ymn,
            "xsiz": xsiz,
            "ysiz": ysiz,
            "zmin": zmin,
            "zmax": zmax,
            "nxdis": nxdis,
            "nydis": nydis,
            "ndmin": ndmin,
            "ndmax": ndmax,
            "radius": radius,
            "nugget": nugget,
            "nst": nst,
            "it1": it1,
            "cc1": cc1,
            "azi1": azi1,
            "hmaj1": hmaj1,
            "hmin1": hmin1,
        }

        self.lateral_sgs_parameters = l_parameters

    def set_g_sgs_GSLIBparameters(self, vertical_simulation_parameters):

        """
        Simulation parameters for growth (vertical) extent simulation

        Parameters
        ----------

        vertical_simulation_parameters: python dictionary with the following parameters -->

        tmin, tmax : all values, regardless of which variable, strictly less than tmin and greater than or equal to tmax are ignored.
        itrans : 0 - no transformation requiered, data already nscored / 1 - for transformation
        ktype : type of interpolation, 0 for simple kriging - 1 for ordinary kriging

        nx, ny : Numbers of blocks. Grid node indices ix, iy increase from 1 to nx, ny respectively (in the positive x,y direction).
        xmn, ymn :  Centre of the first block (xmn,ymn)
        xsiz, ysiz : Size of blocks
        zmin, zmax : The minimum and maximum allowable data values for maximum g simulation (contact opposite to intrusion network).
                    These are used in the back-transformation procedure.
        zmin2, zmax2 : The minimum and maximum allowable data values for minimum g simulation (simultion to condition model to intrusion network).
                    These are used in the back-transformation procedure.
        nxdis, nydis : Number of discretization points for a block. If both nxdis and nydis are set to 1, then point kriging is performed.
        ndmin, ndmax : The minimum and maximum number of original data that should be used to simulate a grid node. If there are fewer than ndmin data points, the node is not simulated.
        radius : The maximum isotropic search radius.

        Returns
        ----
        """
        self._up_to_date = False
        tmin = vertical_simulation_parameters.get("tmin", -9999)
        tmax = vertical_simulation_parameters.get("tmax", 9999)
        itrans = vertical_simulation_parameters.get("itrans", 1)
        ktype = vertical_simulation_parameters.get("ktype", 0)
        nx = vertical_simulation_parameters.get("nx", None)
        ny = vertical_simulation_parameters.get("ny", None)
        xmn = vertical_simulation_parameters.get("xmn", None)
        ymn = vertical_simulation_parameters.get("ymn", None)
        xsiz = vertical_simulation_parameters.get("xsiz", None)
        ysiz = vertical_simulation_parameters.get("ysiz", None)
        zmin = vertical_simulation_parameters.get("zmin", None)
        zmax = vertical_simulation_parameters.get("zmax", None)
        zmin2 = vertical_simulation_parameters.get("zmin2", None)
        zmax2 = vertical_simulation_parameters.get("zmax2", None)
        nxdis = vertical_simulation_parameters.get("nxdis", 1)
        nydis = vertical_simulation_parameters.get("nydis", 1)
        ndmin = vertical_simulation_parameters.get("ndmin", 0)
        ndmax = vertical_simulation_parameters.get("ndmax", 3)
        radius = vertical_simulation_parameters.get("radius", 500)
        nugget = vertical_simulation_parameters.get("nugget", 0)
        nst = vertical_simulation_parameters.get("nst", 1)
        it1 = vertical_simulation_parameters.get("it1", 1)
        cc1 = vertical_simulation_parameters.get("cc1", 1)
        azi1 = vertical_simulation_parameters.get("azi1", 90)
        hmaj1 = vertical_simulation_parameters.get("hmaj1", 999999)
        hmin1 = vertical_simulation_parameters.get("hmin1", 999999)

        g_parameters = {
            "tmin": tmin,
            "tmax": tmax,
            "itrans": itrans,
            "ktype": ktype,
            "nx": nx,
            "ny": ny,
            "xmn": xmn,
            "ymn": ymn,
            "xsiz": xsiz,
            "ysiz": ysiz,
            "zmin": zmin,
            "zmax": zmax,
            "zmin2": zmin2,
            "zmax2": zmax2,
            "nxdis": nxdis,
            "nydis": nydis,
            "ndmin": ndmin,
            "ndmax": ndmax,
            "radius": radius,
            "nugget": nugget,
            "nst": nst,
            "it1": it1,
            "cc1": cc1,
            "azi1": azi1,
            "hmaj1": hmaj1,
            "hmin1": hmin1,
        }

        self.vertical_sgs_parameters = g_parameters

    def make_l_sgs_variogram(self):
        """
        Make variogram for lateral extent simulation
        By default: variogram with no nugget effect, 1 spherical nested structure, isotropic and infinite range

        Parameters
        ----------
        python dictionary with the following parameters -->

        nugget = 0; nst = 1                # nugget effect = 0 and 1 nested structure
        it1 = 1                            # nested structure type (1 - spherical, 2 - exponential, 3 - Gaussian)
        cc1 = 1                            # One structure and no nugget
        azi1 = 90                          # azimuth of this nested structure
        hmaj1 = 9999 ; hmin1 = 9999   # range of this nested structure in the major (hmaj1) and minor (hmin1) direction

        Returns
        ----
        """

        nugget = self.lateral_sgs_parameters.get("nugget")
        nst = self.lateral_sgs_parameters.get("nst")
        it1 = self.lateral_sgs_parameters.get("it1")
        cc1 = self.lateral_sgs_parameters.get("cc1")
        azi1 = self.lateral_sgs_parameters.get("azi1")
        hmaj1 = self.lateral_sgs_parameters.get("hmaj1")
        hmin1 = self.lateral_sgs_parameters.get("hmin1")

        self.lateral_sgs_variogram = GSLIB.make_variogram(
            nugget, nst, it1, cc1, azi1, hmaj1, hmin1
        )

    def make_g_sgs_variogram(self):
        """
        Make variogram for vertical extent simulation
        By default: variogram with no nugget effect, 1 spherical nested structure, isotropic and infinite range

        Parameters
        ----------
        nugget = 0; nst = 1                # nugget effect = 0 and 1 nested structure
        it1 = 1                            # nested structure type (1 - spherical, 2 - exponential, 3 - Gaussian)
        cc1 = 1                            # One structure and no nugget
        azi1 = 90                          # azimuth of this nested structure
        hmaj1 = 9999 ; hmin1 = 9999   # range of this nested structure in the major (hmaj1) and minor (hmin1) direction

        Returns
        ----
        """

        nugget = self.vertical_sgs_parameters.get("nugget")
        nst = self.vertical_sgs_parameters.get("nst")
        it1 = self.vertical_sgs_parameters.get("it1")
        cc1 = self.vertical_sgs_parameters.get("cc1")
        azi1 = self.vertical_sgs_parameters.get("azi1")
        hmaj1 = self.vertical_sgs_parameters.get("hmaj1")
        hmin1 = self.vertical_sgs_parameters.get("hmin1")

        self.vertical_sgs_variogram = GSLIB.make_variogram(
            nugget, nst, it1, cc1, azi1, hmaj1, hmin1
        )

    def set_conceptual_models_parameters(self):

        grid_points_coord1 = self.simulation_grid[2]

        modelcover, minP, maxP, minL, maxL = self.lateral_extent_model()

        if minL == None:
            minL = min(
                self.vertical_contact_data[0]["coord2"].min(),
                self.vertical_contact_data[1]["coord2"].min(),
                self.lateral_contact_data[0]["coord2"].min(),
            )

        if maxL == None:
            maxL = max(
                self.vertical_contact_data[0]["coord2"].max(),
                self.vertical_contact_data[1]["coord2"].max(),
                self.lateral_contact_data[0]["coord2"].max(),
            )

        if minL < 0 and maxL < 0:
            maxL = minL * -1

        if minL > 0 and maxL > 0:
            minL = maxL * -1

        if modelcover == True:
            minP = np.nanmin(grid_points_coord1)
            maxP = np.nanmax(grid_points_coord1)
        else:
            minP = min(
                self.vertical_contact_data[0]["coord1"].min(),
                self.vertical_contact_data[1]["coord1"].min(),
                self.lateral_contact_data[0]["coord1"].min(),
            )
            maxP = max(
                self.vertical_contact_data[0]["coord1"].max(),
                self.vertical_contact_data[1]["coord1"].max(),
                self.lateral_contact_data[0]["coord1"].max(),
            )

        # extra parameters for growth
        mean_growth = self.vertical_contact_data[1].loc[:, "coord0"].mean()
        maxG = self.vertical_contact_data[1]["coord0"].max()
        coord_PL_for_maxG = (
            self.vertical_contact_data[1][
                self.vertical_contact_data[1].coord0
                == self.vertical_contact_data[1].coord0.max()
            ]
            .loc[:, ["coord1", "coord2"]]
            .to_numpy()
        )

        vertex = [coord_PL_for_maxG[0][0], coord_PL_for_maxG[0][1], maxG]

        # conceptual_model_parameters = {}
        self.conceptual_model_parameters["minP"] = minP
        self.conceptual_model_parameters["maxP"] = maxP
        self.conceptual_model_parameters["minL"] = minL
        self.conceptual_model_parameters["maxL"] = maxL
        self.conceptual_model_parameters["model_cover"] = modelcover
        self.conceptual_model_parameters["mean_growth"] = mean_growth
        self.conceptual_model_parameters["vertex"] = vertex

    def simulate_lateral_thresholds_original(self):
        """
        Simulate residual values along L frame axis,
        and compute l thresholds for each point of a pre-defined grid

        Parameters
        ----------

        Returns
        -------
        """

        # -- get grid points and evaluated values in the intrusion frame
        grid_points_coord1 = self.simulation_grid[2]

        self.set_conceptual_models_parameters()

        minP = self.conceptual_model_parameters.get("minP")
        maxP = self.conceptual_model_parameters.get("maxP")
        minL = self.conceptual_model_parameters.get("minL")
        maxL = self.conceptual_model_parameters.get("maxL")

        if self.width_data[0] == False:  # i.e., no lateral data for side L<0
            print(
                "Not enought lateral data for simulation of side L<0, Using roof/floor data to condition the conceptual model"
            )

            # -- try using vertical data to set some points and run SGS
            vertical_data = pd.concat(
                [self.vertical_contact_data[0], self.vertical_contact_data[1]]
            )
            vertical_data.loc[
                :, ["conceptual_maxside", "conceptual_minside"]
            ] = self.lateral_extent_model(
                lateral_contact_data=vertical_data,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )

            data_minL_temp = vertical_data[vertical_data["coord2"] < 0].copy()

            inputsimdata_minL = (
                data_minL_temp[
                    data_minL_temp["coord2"] <= data_minL_temp["conceptual_minside"]
                ]
                .loc[
                    :,
                    ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_minside"],
                ]
                .copy()
            )
            inputsimdata_minL.loc[:, "l_residual"] = (
                inputsimdata_minL.loc[:, "conceptual_minside"]
                - inputsimdata_minL.loc[:, "coord2"]
            )
            inputsimdata_minL.rename(
                columns={"conceptual_minside": "l_conceptual"}, inplace=True
            )
            inputsimdata_minL.reset_index(inplace=True)

            if len(inputsimdata_minL) < 3:
                # create random points along coordinate 1, and evaluate them in conceptual model. Residual = 0. Add points to dataframe containing input for sgs
                print(
                    "Simulation of lateral side L<0: No enought roof/floor data to condition the conceptual model, lateral contact equivalent to conceptual"
                )

                random_p = pd.DataFrame(
                    np.random.randint(minP, maxP, 10), columns=["coord1"]
                )
                conceptual_l = self.lateral_extent_model(
                    lateral_contact_data=random_p,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                inputsimdata_minL_ = pd.DataFrame(
                    np.vstack(
                        [conceptual_l[:, 1], random_p.loc[:, "coord1"].to_numpy()]
                    ).T,
                    columns=["l_conceptual", "coord1"],
                )
                inputsimdata_minL_.loc[:, "l_residual"] = 0
                inputsimdata_minL = pd.concat([inputsimdata_minL, inputsimdata_minL_])

            inputsimdata_minL.loc[:, "ref_coord"] = 0

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
            inputsimdata_minL = data_minL.loc[
                :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
            ].copy()
            inputsimdata_minL.loc[:, "l_residual"] = data_residual_minL
            inputsimdata_minL.loc[:, "l_conceptual"] = data_conceptual_minL[:, 1]

            inputsimdata_minL.reset_index(inplace=True)
            inputsimdata_minL.loc[:, "ref_coord"] = 0

        if self.width_data[1] == False:  # i.e., no lateral data for side L>0
            print(
                "Not enought lateral data for simulation of side L>=0, Using roof/floor data to condition the conceptual model"
            )

            # -- try using vertical data to set some points and run SGS
            vertical_data = pd.concat(
                [self.vertical_contact_data[0], self.vertical_contact_data[1]]
            )
            vertical_data.loc[
                :, ["conceptual_maxside", "conceptual_minside"]
            ] = self.lateral_extent_model(
                lateral_contact_data=vertical_data,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )

            data_maxL_temp = vertical_data[vertical_data["coord2"] >= 0].copy()

            inputsimdata_maxL = (
                data_maxL_temp[
                    data_maxL_temp["coord2"] >= data_maxL_temp["conceptual_maxside"]
                ]
                .loc[
                    :,
                    ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_maxside"],
                ]
                .copy()
            )
            inputsimdata_maxL.loc[:, "l_residual"] = (
                inputsimdata_maxL.loc[:, "conceptual_maxside"]
                - inputsimdata_maxL.loc[:, "coord2"]
            )
            inputsimdata_maxL.rename(
                columns={"conceptual_maxside": "l_conceptual"}, inplace=True
            )
            inputsimdata_maxL.reset_index(inplace=True)

            if len(inputsimdata_maxL) < 3:
                print(
                    "Simulation of lateral side L>=0: No enought roof/floor data to condition the conceptual model, lateral contact equivalent to conceptual"
                )
                # create random points along coordinate 1, and evaluate them conceptual model. Residual = 0. Add points to dataframe containing input for sgs
                random_p = pd.DataFrame(
                    np.random.randint(minP, maxP, 10), columns=["coord1"]
                )
                conceptual_l = self.lateral_extent_model(
                    lateral_contact_data=random_p,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                inputsimdata_maxL_ = pd.DataFrame(
                    np.vstack(
                        [conceptual_l[:, 0], random_p.loc[:, "coord1"].to_numpy()]
                    ).T,
                    columns=["l_conceptual", "coord1"],
                )
                inputsimdata_maxL_.loc[:, "l_residual"] = 0
                inputsimdata_maxL = pd.concat([inputsimdata_maxL, inputsimdata_maxL_])

            inputsimdata_maxL.loc[:, "ref_coord"] = 0

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
            inputsimdata_maxL = data_maxL.loc[
                :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
            ].copy()
            inputsimdata_maxL.loc[:, "l_residual"] = data_residual_maxL
            inputsimdata_maxL.loc[:, "l_conceptual"] = data_conceptual_maxL[:, 0]

            inputsimdata_maxL.reset_index(inplace=True)
            inputsimdata_maxL.loc[:, "ref_coord"] = 0

        self.lateral_sgs_input_data = [inputsimdata_minL, inputsimdata_maxL]

        # -- Compute simulation parameters if not defined
        # ---- compute lower and uper fence of simulated values using quartiles

        p25_lmin = np.percentile(
            inputsimdata_minL.copy().drop_duplicates(subset="l_residual").l_residual, 25
        )
        p75_lmin = np.percentile(
            inputsimdata_minL.copy().drop_duplicates(subset="l_residual").l_residual, 75
        )
        p25_lmax = np.percentile(
            inputsimdata_maxL.copy().drop_duplicates(subset="l_residual").l_residual, 25
        )
        p75_lmax = np.percentile(
            inputsimdata_maxL.copy().drop_duplicates(subset="l_residual").l_residual, 75
        )

        if self.lateral_sgs_parameters.get("zmin") == None:
            self.lateral_sgs_parameters["zmin_lmin"] = p25_lmin - 1.5 * (
                p75_lmin - p25_lmin
            )
            self.lateral_sgs_parameters["zmin_lmax"] = p25_lmax - 1.5 * (
                p75_lmax - p25_lmax
            )
        else:
            self.lateral_sgs_parameters["zmin_lmin"] = self.lateral_sgs_parameters.get(
                "zmin"
            )
            self.lateral_sgs_parameters["zmin_lmax"] = self.lateral_sgs_parameters.get(
                "zmin"
            )

        if self.lateral_sgs_parameters.get("zmax") == None:
            self.lateral_sgs_parameters["zmax_lmin"] = p75_lmin + 1.5 * (
                p75_lmin - p25_lmin
            )
            self.lateral_sgs_parameters["zmax_lmax"] = p75_lmax + 1.5 * (
                p75_lmax - p25_lmax
            )

        else:
            self.lateral_sgs_parameters["zmax_lmin"] = self.lateral_sgs_parameters.get(
                "zmax"
            )
            self.lateral_sgs_parameters["zmax_lmax"] = self.lateral_sgs_parameters.get(
                "zmax"
            )

        if self.lateral_sgs_parameters.get("xmn") == None:
            self.lateral_sgs_parameters["xmn"] = np.nanmin(grid_points_coord1)
        if self.lateral_sgs_parameters.get("xsiz") == None:
            nx = self.lateral_sgs_parameters.get("nx")
            minC0 = np.nanmin(grid_points_coord1)
            maxC0 = np.nanmax(grid_points_coord1)
            self.lateral_sgs_parameters["xsiz"] = (maxC0 - minC0) / nx

        # -- Simulation of lateral extent
        tmin = self.lateral_sgs_parameters.get("tmin")
        tmax = self.lateral_sgs_parameters.get("tmax")
        itrans = self.lateral_sgs_parameters.get("itrans")
        ktype = self.lateral_sgs_parameters.get("ktype")
        nx = self.lateral_sgs_parameters.get("nx")
        ny = self.lateral_sgs_parameters.get("ny")
        xmn = self.lateral_sgs_parameters.get("xmn")
        ymn = self.lateral_sgs_parameters.get("ymn")
        xsiz = self.lateral_sgs_parameters.get("xsiz")
        ysiz = self.lateral_sgs_parameters.get("ysiz")
        zmin = self.lateral_sgs_parameters.get("zmin_lmin")
        zmax = self.lateral_sgs_parameters.get("zmax_lmin")
        nxdis = self.lateral_sgs_parameters.get("nxdis")
        nydis = self.lateral_sgs_parameters.get("nydis")
        ndmin = self.lateral_sgs_parameters.get("ndmin")
        ndmax = self.lateral_sgs_parameters.get("ndmax")
        radius = self.lateral_sgs_parameters.get("radius")

        # self.lateral_sgs_input_data = [inputsimdata_minL, inputsimdata_maxL]

        l_min_simulation = geostats.sgsim(
            self.lateral_sgs_input_data[0],
            "coord1",
            "ref_coord",
            "l_residual",
            wcol=-1,
            scol=-1,
            tmin=tmin,
            tmax=tmax,
            itrans=itrans,
            ismooth=0,
            dftrans=0,
            tcol=0,
            twtcol=0,
            zmin=zmin,
            zmax=zmax,
            ltail=1,
            ltpar=0.0,
            utail=1,
            utpar=0.3,
            nsim=1,
            nx=nx,
            xmn=xmn,
            xsiz=xsiz,
            ny=ny,
            ymn=ymn,
            ysiz=ysiz,
            seed=73073,
            ndmin=ndmin,
            ndmax=ndmax,
            nodmax=1,
            mults=0,
            nmult=2,
            noct=-1,
            radius=radius,
            radius1=10,
            sang1=0,
            mxctx=1,
            mxcty=1,
            ktype=ktype,
            colocorr=0.0,
            sec_map=0,
            vario=self.lateral_sgs_variogram,
        )

        zmin = self.lateral_sgs_parameters.get("zmin_lmax")
        zmax = self.lateral_sgs_parameters.get("zmax_lmax")

        zmin = self.lateral_sgs_parameters.get("zmin_lmax")
        zmax = self.lateral_sgs_parameters.get("zmax_lmax")

        l_max_simulation = geostats.sgsim(
            self.lateral_sgs_input_data[1],
            "coord1",
            "ref_coord",
            "l_residual",
            wcol=-1,
            scol=-1,
            tmin=tmin,
            tmax=tmax,
            itrans=itrans,
            ismooth=0,
            dftrans=0,
            tcol=0,
            twtcol=0,
            zmin=zmin,
            zmax=zmax,
            ltail=1,
            ltpar=0.0,
            utail=1,
            utpar=0.3,
            nsim=1,
            nx=nx,
            xmn=xmn,
            xsiz=xsiz,
            ny=ny,
            ymn=ymn,
            ysiz=ysiz,
            seed=73073,
            ndmin=ndmin,
            ndmax=ndmax,
            nodmax=1,
            mults=0,
            nmult=2,
            noct=-1,
            radius=radius,
            radius1=10,
            sang1=0,
            mxctx=1,
            mxcty=1,
            ktype=ktype,
            colocorr=0.0,
            sec_map=0,
            vario=self.lateral_sgs_variogram,
        )

        self.simulationGSLIB_s_outcome = [l_min_simulation, l_max_simulation]

        # -- Create dataframe containing S threshold for each grid point

        propagation_grid_model = np.linspace(xmn, xmn + (nx * xsiz), nx)

        lateral_thresholds = pd.DataFrame(
            columns=[
                "coord1",
                "min_l_residual",
                "min_l_threshold",
                "max_l_residual",
                "max_l_threshold",
            ]
        )

        lateral_thresholds["coord1"] = propagation_grid_model

        model_conceptual_l = self.lateral_extent_model(
            lateral_contact_data=lateral_thresholds,
            minP=minP,
            maxP=maxP,
            minS=minL,
            maxS=maxL,
        )
        lateral_thresholds["conceptual_minl"] = model_conceptual_l[:, 1]
        lateral_thresholds["min_l_residual"] = l_min_simulation[1]
        lateral_thresholds["min_l_threshold"] = (
            model_conceptual_l[:, 1] - l_min_simulation[1]
        )

        lateral_thresholds["conceptual_maxl"] = model_conceptual_l[:, 0]
        lateral_thresholds["max_l_residual"] = l_max_simulation[1]
        lateral_thresholds["max_l_threshold"] = (
            model_conceptual_l[:, 0] - l_max_simulation[1]
        )
        lateral_thresholds.sort_values(["coord1"], ascending=[True], inplace=True)

        # Ignore simulated data outside area covered by input data

        lateral_thresholds.loc[
            lateral_thresholds.coord1 < minP, ["min_l_threshold", "max_l_threshold"]
        ] = [0.00001, 0.00001]

        lateral_thresholds.loc[
            lateral_thresholds.coord1 > maxP, ["min_l_threshold", "max_l_threshold"]
        ] = [0.00001, 0.00001]

        self.lateral_simulated_thresholds = lateral_thresholds

    def simulate_lateral_thresholds(self):
        """
        Simulate residual values along L frame axis,
        and compute l thresholds for each point of a pre-defined grid

        Parameters
        ----------

        Returns
        -------
        """

        # -- get grid points and evaluated values in the intrusion frame

        # grid_points_coord1 = self.simulation_grid[2]

        # -- generate data frame containing input data for simulation

        self.set_conceptual_models_parameters()

        minP = self.conceptual_model_parameters.get("minP")
        maxP = self.conceptual_model_parameters.get("maxP")
        minL = self.conceptual_model_parameters.get("minL")
        maxL = self.conceptual_model_parameters.get("maxL")

        if self.width_data[0] == False:  # i.e., no lateral data for side L<0
            print(
                "Not enought lateral data for simulation of side L<0, Using roof/floor data to condition the conceptual model"
            )

            # -- try using vertical data to set some points and run SGS
            vertical_data = pd.concat(
                [self.vertical_contact_data[0], self.vertical_contact_data[1]]
            )
            vertical_data.loc[
                :, ["conceptual_maxside", "conceptual_minside"]
            ] = self.lateral_extent_model(
                lateral_contact_data=vertical_data,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )

            data_minL_temp = vertical_data[vertical_data["coord2"] < 0].copy()

            inputsimdata_minL = (
                data_minL_temp[
                    data_minL_temp["coord2"] <= data_minL_temp["conceptual_minside"]
                ]
                .loc[
                    :,
                    ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_minside"],
                ]
                .copy()
            )
            inputsimdata_minL.loc[:, "l_residual"] = (
                inputsimdata_minL.loc[:, "conceptual_minside"]
                - inputsimdata_minL.loc[:, "coord2"]
            )
            inputsimdata_minL.rename(
                columns={"conceptual_minside": "l_conceptual"}, inplace=True
            )
            inputsimdata_minL.reset_index(inplace=True)

            if len(inputsimdata_minL) < 3:
                # create random points along coordinate 1, and evaluate them in conceptual model. Residual = 0. Add points to dataframe containing input for sgs
                print(
                    "Simulation of lateral side L<0: No enought roof/floor data to condition the conceptual model, lateral contact equivalent to conceptual"
                )

                random_p = pd.DataFrame(
                    np.random.randint(minP, maxP, 10), columns=["coord1"]
                )
                conceptual_l = self.lateral_extent_model(
                    lateral_contact_data=random_p,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                inputsimdata_minL_ = pd.DataFrame(
                    np.vstack(
                        [conceptual_l[:, 1], random_p.loc[:, "coord1"].to_numpy()]
                    ).T,
                    columns=["l_conceptual", "coord1"],
                )
                inputsimdata_minL_.loc[:, "l_residual"] = 0
                inputsimdata_minL = pd.concat([inputsimdata_minL, inputsimdata_minL_])

            inputsimdata_minL.loc[:, "ref_coord"] = 0

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
            inputsimdata_minL = data_minL.loc[
                :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
            ].copy()
            inputsimdata_minL.loc[:, "l_residual"] = data_residual_minL
            inputsimdata_minL.loc[:, "l_conceptual"] = data_conceptual_minL[:, 1]

            inputsimdata_minL.reset_index(inplace=True)
            inputsimdata_minL.loc[:, "ref_coord"] = 0

        if self.width_data[1] == False:  # i.e., no lateral data for side L>0
            print(
                "Not enought lateral data for simulation of side L>=0, Using roof/floor data to condition the conceptual model"
            )

            # -- try using vertical data to set some points and run SGS
            vertical_data = pd.concat(
                [self.vertical_contact_data[0], self.vertical_contact_data[1]]
            )
            vertical_data.loc[
                :, ["conceptual_maxside", "conceptual_minside"]
            ] = self.lateral_extent_model(
                lateral_contact_data=vertical_data,
                minP=minP,
                maxP=maxP,
                minS=minL,
                maxS=maxL,
            )

            data_maxL_temp = vertical_data[vertical_data["coord2"] >= 0].copy()

            inputsimdata_maxL = (
                data_maxL_temp[
                    data_maxL_temp["coord2"] >= data_maxL_temp["conceptual_maxside"]
                ]
                .loc[
                    :,
                    ["X", "Y", "Z", "coord0", "coord1", "coord2", "conceptual_maxside"],
                ]
                .copy()
            )
            inputsimdata_maxL.loc[:, "l_residual"] = (
                inputsimdata_maxL.loc[:, "conceptual_maxside"]
                - inputsimdata_maxL.loc[:, "coord2"]
            )
            inputsimdata_maxL.rename(
                columns={"conceptual_maxside": "l_conceptual"}, inplace=True
            )
            inputsimdata_maxL.reset_index(inplace=True)

            if len(inputsimdata_maxL) < 3:
                print(
                    "Simulation of lateral side L>=0: No enought roof/floor data to condition the conceptual model, lateral contact equivalent to conceptual"
                )
                # create random points along coordinate 1, and evaluate them conceptual model. Residual = 0. Add points to dataframe containing input for sgs
                random_p = pd.DataFrame(
                    np.random.randint(minP, maxP, 10), columns=["coord1"]
                )
                conceptual_l = self.lateral_extent_model(
                    lateral_contact_data=random_p,
                    minP=minP,
                    maxP=maxP,
                    minS=minL,
                    maxS=maxL,
                )
                inputsimdata_maxL_ = pd.DataFrame(
                    np.vstack(
                        [conceptual_l[:, 0], random_p.loc[:, "coord1"].to_numpy()]
                    ).T,
                    columns=["l_conceptual", "coord1"],
                )
                inputsimdata_maxL_.loc[:, "l_residual"] = 0
                inputsimdata_maxL = pd.concat([inputsimdata_maxL, inputsimdata_maxL_])

            inputsimdata_maxL.loc[:, "ref_coord"] = 0

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
            inputsimdata_maxL = data_maxL.loc[
                :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
            ].copy()
            inputsimdata_maxL.loc[:, "l_residual"] = data_residual_maxL
            inputsimdata_maxL.loc[:, "l_conceptual"] = data_conceptual_maxL[:, 0]

            inputsimdata_maxL.reset_index(inplace=True)
            inputsimdata_maxL.loc[:, "ref_coord"] = 0

        self.lateral_sgs_input_data = [inputsimdata_minL, inputsimdata_maxL]

    def interpolate_lateral_thresholds(self):
        return "TO DO"

    def simulate_growth_thresholds_original(self):
        """
        Simulate g residual values and compute g thresholds for each point of a pre-defined grid.
        Computes two sets of g thresholds: one to constraint the contact opposite to the intrusion network (using residual values),
        and another one to better condition the intrusion network contact to the data.

        Parameters
        ----------

        Returns
        -------
        """

        # -- Get grid points and evaluated values in the intrusion frame
        grid_points = self.simulation_grid
        grid_points_coord0 = self.simulation_grid[1]
        grid_points_coord1 = self.simulation_grid[2]
        grid_points_coord2 = self.simulation_grid[3]

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
        data_residual_G = (
            data_conceptual_G[:, 1] - other_contact_data.loc[:, "coord0"]
        ).to_numpy()
        inputsimdata_maxG = other_contact_data.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()
        inputsimdata_maxG.loc[:, "g_residual"] = data_residual_G
        inputsimdata_maxG.loc[:, "g_conceptual"] = data_conceptual_G[:, 1]

        # --- growth simulation input data for intrusion network conditioning
        inputsimdata_inetG = inet_data.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()

        self.vertical_sgs_input_data = [inputsimdata_maxG, inputsimdata_inetG]

        # Simulation
        # --- compute simulation parameters if not defined
        if self.vertical_sgs_parameters.get("nx") == None:
            self.vertical_sgs_parameters["nx"] = grid_points[4][0]
        if self.vertical_sgs_parameters.get("ny") == None:
            self.vertical_sgs_parameters["ny"] = grid_points[4][1]
        if self.vertical_sgs_parameters.get("xmn") == None:
            self.vertical_sgs_parameters["xmn"] = np.nanmin(grid_points_coord1)
        if self.vertical_sgs_parameters.get("ymn") == None:
            self.vertical_sgs_parameters["ymn"] = np.nanmin(grid_points_coord2)
        if self.vertical_sgs_parameters.get("xsiz") == None:
            nx = self.vertical_sgs_parameters.get("nx")
            minC1 = np.nanmin(grid_points_coord1)
            maxC1 = np.nanmax(grid_points_coord1)
            self.vertical_sgs_parameters["xsiz"] = (maxC1 - minC1) / nx
        if self.vertical_sgs_parameters.get("ysiz") == None:
            ny = self.vertical_sgs_parameters.get("ny")
            minC2 = np.nanmin(grid_points_coord2)
            maxC2 = np.nanmax(grid_points_coord2)
            self.vertical_sgs_parameters["ysiz"] = (maxC2 - minC2) / ny
        if self.vertical_sgs_parameters.get("zmin") == None:
            self.vertical_sgs_parameters["zmin"] = inputsimdata_maxG["g_residual"].min()
        if self.vertical_sgs_parameters.get("zmax") == None:
            self.vertical_sgs_parameters["zmax"] = inputsimdata_maxG["g_residual"].max()
        if self.vertical_sgs_parameters.get("zmin2") == None:
            self.vertical_sgs_parameters["zmin2"] = inputsimdata_inetG["coord0"].min()
        if self.vertical_sgs_parameters.get("zmax2") == None:
            self.vertical_sgs_parameters["zmax2"] = inputsimdata_inetG["coord0"].max()

        # --- get parameters for simulation
        tmin = self.vertical_sgs_parameters.get("tmin")
        tmax = self.vertical_sgs_parameters.get("tmax")
        itrans = self.vertical_sgs_parameters.get("itrans")
        ktype = self.vertical_sgs_parameters.get("ktype")
        nx = self.vertical_sgs_parameters.get("nx")
        ny = self.vertical_sgs_parameters.get("ny")
        xmn = self.vertical_sgs_parameters.get("xmn")
        ymn = self.vertical_sgs_parameters.get("ymn")
        xsiz = self.vertical_sgs_parameters.get("xsiz")
        ysiz = self.vertical_sgs_parameters.get("ysiz")
        zmin = self.vertical_sgs_parameters.get("zmin")
        zmax = self.vertical_sgs_parameters.get("zmax")
        zmin2 = self.vertical_sgs_parameters.get("zmin2")
        zmax2 = self.vertical_sgs_parameters.get("zmax2")
        nxdis = self.vertical_sgs_parameters.get("nxdis")
        nydis = self.vertical_sgs_parameters.get("nydis")
        ndmin = self.vertical_sgs_parameters.get("ndmin")
        ndmax = self.vertical_sgs_parameters.get("ndmax")
        radius = self.vertical_sgs_parameters.get("radius")

        # --- simulation of residual values related to other contact( opposite to intrusion network)
        g_max_simulation = geostats.sgsim(
            inputsimdata_maxG,
            "coord1",
            "coord2",
            "g_residual",
            wcol=-1,
            scol=-1,
            tmin=tmin,
            tmax=tmax,
            itrans=itrans,
            ismooth=0,
            dftrans=0,
            tcol=0,
            twtcol=0,
            zmin=zmin,
            zmax=zmax,
            ltail=1,
            ltpar=0.0,
            utail=1,
            utpar=0.3,
            nsim=1,
            nx=nx,
            xmn=xmn,
            xsiz=xsiz,
            ny=ny,
            ymn=ymn,
            ysiz=ysiz,
            seed=73073,
            ndmin=ndmin,
            ndmax=ndmax,
            nodmax=1,
            mults=0,
            nmult=2,
            noct=-1,
            radius=radius,
            radius1=10,
            sang1=0,
            mxctx=1,
            mxcty=1,
            ktype=ktype,
            colocorr=0.0,
            sec_map=0,
            vario=self.vertical_sgs_variogram,
        )

        # --- simulation to improve conditioning of intrusion network contact
        g_min_simulation = geostats.sgsim(
            inputsimdata_inetG,
            "coord1",
            "coord2",
            "coord0",
            wcol=-1,
            scol=-1,
            tmin=tmin,
            tmax=tmax,
            itrans=itrans,
            ismooth=0,
            dftrans=0,
            tcol=0,
            twtcol=0,
            zmin=zmin2,
            zmax=zmax2,
            ltail=1,
            ltpar=0.0,
            utail=1,
            utpar=0.3,
            nsim=1,
            nx=nx,
            xmn=xmn,
            xsiz=xsiz,
            ny=ny,
            ymn=ymn,
            ysiz=ysiz,
            seed=73073,
            ndmin=ndmin,
            ndmax=ndmax,
            nodmax=1,
            mults=0,
            nmult=2,
            noct=-1,
            radius=radius,
            radius1=10,
            sang1=0,
            mxctx=1,
            mxcty=1,
            ktype=ktype,
            colocorr=0.0,
            sec_map=0,
            vario=self.vertical_sgs_variogram,
        )

        # self.simulationGSLIB_g_outcome = [g_min_simulation, g_max_simulation]
        # -- Create dataframe containing S threshold for each grid point
        lower_extent_gpl = [0, xmn, ymn]
        upper_extent_gpl = [0, xmn + xsiz * nx, ymn + ysiz * ny]

        # -- transform SGS output array to Intrusion Frame coordinates
        # -- grid_from_array returns a matrix of [node_i,node_j,g,p,s,values in array]
        g_residual = grid_from_array(
            g_max_simulation, ["X", 0], lower_extent_gpl, upper_extent_gpl
        )
        g_minimum = grid_from_array(
            g_min_simulation, ["X", 0], lower_extent_gpl, upper_extent_gpl
        )

        simulation_g_threshold = pd.DataFrame(
            columns=[
                "coord0",
                "coord1",
                "coord2",
                "g_residual",
                "g_maximum",
                "g_minimum",
            ]
        )
        simulation_g_threshold.loc[:, "coord0"] = g_residual[:, 2]
        simulation_g_threshold.loc[:, "coord1"] = g_residual[:, 3]
        simulation_g_threshold.loc[:, "coord2"] = g_residual[:, 4]
        simulation_g_threshold.loc[:, "g_residual"] = g_residual[:, 5]
        simulation_g_threshold.loc[:, "g_minimum"] = g_minimum[:, 5]
        conceptual_maxg = self.vertical_extent_model(
            simulation_g_threshold,
            mean_growth=meanG,
            minP=minP,
            maxP=maxP,
            minS=minL,
            maxS=maxL,
            vertex=vertex,
        )
        g_maximum = conceptual_maxg[:, 1] - g_residual[:, 5]
        simulation_g_threshold.loc[:, "g_maximum"] = g_maximum

        grid_for_growth_evaluation = StructuredGrid2D(
            origin=[xmn, ymn], nsteps=[nx, ny], step_vector=[xsiz, ysiz]
        )
        grid_for_growth_evaluation.update_property("g_maximum", g_maximum)
        grid_for_growth_evaluation.update_property("g_minimum", g_minimum[:, 5])

        self.growth_simulated_thresholds = simulation_g_threshold
        self.growth_simulated_thresholds_grid = grid_for_growth_evaluation

    def simulate_growth_thresholds(self):
        """
        Simulate g residual values and compute g thresholds for each point of a pre-defined grid.
        Computes two sets of g thresholds: one to constraint the contact opposite to the intrusion network (using residual values),
        and another one to better condition the intrusion network contact to the data.

        Parameters
        ----------

        Returns
        -------
        """

        # -- Get grid points and evaluated values in the intrusion frame
        grid_points = self.simulation_grid
        grid_points_coord0 = self.simulation_grid[1]
        grid_points_coord1 = self.simulation_grid[2]
        grid_points_coord2 = self.simulation_grid[3]

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
        data_residual_G = (
            data_conceptual_G[:, 1] - other_contact_data.loc[:, "coord0"]
        ).to_numpy()
        inputsimdata_maxG = other_contact_data.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()
        inputsimdata_maxG.loc[:, "g_residual"] = data_residual_G
        inputsimdata_maxG.loc[:, "g_conceptual"] = data_conceptual_G[:, 1]

        # --- growth simulation input data for intrusion network conditioning
        inputsimdata_inetG = inet_data.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()

        self.vertical_sgs_input_data = [inputsimdata_maxG, inputsimdata_inetG]

        # self.growth_simulated_thresholds = simulation_g_threshold
        # self.growth_simulated_thresholds_grid = grid_for_growth_evaluation

    def build(
        self,
        vertical_extent_sgs_parameters={},
        lateral_extent_sgs_parameters={},
        geometric_scaling_parameters={},
        **kwargs,
    ):
        """Main building function for intrusion. Calculates variogram and simulate thresholds along frame axes

        Parameters
        ----------
        vertical_extent_sgs_parameters : dict, optional
            parameters for the vertical sequential gaussian simulation, by default {}
        lateral_extent_sgs_parameters : dict, optional
            parameters for the vertical sequential gaussian simulation, by default {}
        """
        self.prepare_data(geometric_scaling_parameters)
        self.create_grid_for_simulation()
        self.set_l_sgs_GSLIBparameters(lateral_extent_sgs_parameters)
        self.set_g_sgs_GSLIBparameters(vertical_extent_sgs_parameters)
        self.make_l_sgs_variogram()
        self.make_g_sgs_variogram()
        self.simulate_lateral_thresholds()
        self.simulate_growth_thresholds()
        self.feature.growth_simulated_thresholds = (
            []
        )  # self.growth_simulated_thresholds
        self.feature.lateral_simulated_thresholds = (
            []
        )  # self.lateral_simulated_thresholds
        self.feature.growth_simulated_thresholds_grid = (
            self.growth_simulated_thresholds_grid
        )

    def update(self):
        self.build(**self.build_arguments)
        self._up_to_date = True

    def up_to_date(self, callback=None):
        """
        check if the feature is uptodate
        if its not update.

        Parameters
        ----------
        callback : function
            a function that is called when the feature is updated

        """

        # for f in self.ts:
        #     f.builder.up_to_date(callback=callback)
        # has anything changed in the builder since we built the feature? if so update
        if self._up_to_date == False:
            self.update()
            if callable(callback):
                callback(1)
            return
        if callable(callback):
            callback(1)
