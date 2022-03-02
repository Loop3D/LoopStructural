import numpy as np
import pandas as pd

# import logging
from LoopStructural.utils import getLogger
logger = getLogger(__name__)

# import GSLIB library
try:
    import geostatspy.GSLIB as GSLIB  # GSLIB utilities, viz and wrapped functions
    import geostatspy.geostats as geostats  # GSLIB converted to Python

except ImportError:
    logger.error("GeostatPy not installed \n" "pip install geostatspy")

from LoopStructural.modelling.intrusions.intrusion_support_functions import *


class IntrusionFeature:
    """
    Intrusion feature is a class to represent an intrusion, using a distance scalar field to the intrusion contact.
    """

    def __init__(self, name, model=None):

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
        self.model = model
        self.intrusion_frame = None
        self.builder = None

        # conceptual models
        self.lateral_extent_model = None
        self.vertical_extent_model = None

        # contact data
        self.lateral_contact_data = None
        self.vertical_contact_data = None

        # sequential gaussian simulation parameters
        self.lateral_sgs_parameters = None
        self.vertical_sgs_parameters = None

        self.lateral_sgs_variogram = None
        self.vertical_sgs_variogram = None

        self.lateral_sgs_input_data = None
        self.vertical_sgs_input_data = None

        # simulated thresholds
        self.lateral_simulated_thresholds = None

        # self.simulation_growth_data = None
        # self.simulation_lateral_data = None
        # self.inet_poly = inet_polygon

        # self.intrusion_feature_network = None
        # self.intrusion_feature_frame = None
        # self.intrusion_feature_body = None
        # self.intrusion_builder = None

        # self.intrusion_indicator_function = None
        # self.evaluated_points = None

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

    def set_intrusion_builder(self, builder):
        self.intrusion_builder = builder

    def set_intrusion_network(self, intrusion_network):
        self.intrusion_feature_network = intrusion_network

    def set_intrusion_frame(self, intrusion_frame):
        self.intrusion_feature_frame = intrusion_frame

    def set_intrusion_body(self, intrusion_body):
        self.intrusion_feature_body = intrusion_body

    def set_data_for_extent_simulation(self, intrusion_data):
        """Set data for lateral extent (distances in L axis) simulation.
        Data points are evaluated in the intrusion frame coordinates,
        and separated between each sides (> or <0). Evaluated data points are then used to simulate thresholds distances

        Set data for vertical extent (distances in G axis) simulation
        Data points are evaluated in the intrusion frame coordinates,
        and separated between roof and floor contact. 
        Evaluated data points are then used to simulate thresholds distances.


        Parameters
        ----------

        """

        data = intrusion_data.copy()
        data_xyz = data.loc[:, ["X", "Y", "Z"]].to_numpy()
        data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(data_xyz)
        data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(data_xyz)
        data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(data_xyz)
        data_minside = data[
            (data["intrusion_side"] == True) & (data["coord2"] <= 0)
        ].copy()
        data_minside.reset_index(inplace=True)

        data_maxside = data[
            (data["intrusion_side"] == True) & (data["coord2"] > 0)
        ].copy()
        data_maxside.reset_index(inplace=True)

        data_sides = pd.concat([data_minside, data_maxside])
        data_sides.reset_index(inplace=True)

        self.lateral_contact_data = [data_sides, data_minside, data_maxside]
        # return data_sides, data_minside, data_maxside

        intrusion_network_data_xyz = self.builder.intrusion_network_data.loc[:, ["X", "Y", "Z"]].to_numpy()
        intrusion_network_data = self.builder.intrusion_network_data.loc[:, ["X", "Y", "Z"]].copy()
        intrusion_network_data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.reset_index(inplace=True)

        other_contact_data_xyz = self.builder.other_contact_data.loc[:, ["X", "Y", "Z"]].to_numpy()
        other_contact_data = self.builder.other_contact_data.loc[:, ["X", "Y", "Z"]].copy()
        other_contact_data.loc[:, "coord0"] = self.intrusion_frame[0].evaluate_value(other_contact_data_xyz)
        other_contact_data.loc[:, "coord1"] = self.intrusion_frame[1].evaluate_value(other_contact_data_xyz)
        other_contact_data.loc[:, "coord2"] = self.intrusion_frame[2].evaluate_value(other_contact_data_xyz)
        other_contact_data.reset_index(inplace=True)

        self.vertical_contact_data = [intrusion_network_data, other_contact_data]

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
        s_parameters : dictionary of parameters to be used for lateral extent simulation. Parameters definition above
        ----
        """

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
        }

        self.vertical_sgs_parameters = g_parameters
    
    def make_l_sgs_variogram(self, lateral_simulation_parameters):
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
        nugget = lateral_simulation_parameters.get("nugget", 0)
        nst = lateral_simulation_parameters.get("nst", 1)
        it1 = lateral_simulation_parameters.get("it1", 1)
        cc1 = lateral_simulation_parameters.get("cc1", 1)
        azi1 = lateral_simulation_parameters.get("azi1", 90)
        hmaj1 = lateral_simulation_parameters.get("hmaj1", 999999)
        hmin1 = lateral_simulation_parameters.get("hmin1", 999999)

        variogram = GSLIB.make_variogram(nugget, nst, it1, cc1, azi1, hmaj1, hmin1)
        self.lateral_sgs_variogram = variogram

    def make_g_sgs_variogram(self, vertical_simulation_parameters):
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

        nugget = vertical_simulation_parameters.get("nugget", 0)
        nst = vertical_simulation_parameters.get("nst", 1)
        it1 = vertical_simulation_parameters.get("it1", 1)
        cc1 = vertical_simulation_parameters.get("cc1", 1)
        azi1 = vertical_simulation_parameters.get("azi1", 90)
        hmaj1 = vertical_simulation_parameters.get("hmaj1", 999999)
        hmin1 = vertical_simulation_parameters.get("hmin1", 999999)

        variogram = GSLIB.make_variogram(nugget, nst, it1, cc1, azi1, hmaj1, hmin1)
        self.vertical_sgs_variogram = variogram

    def simulate_lateral_thresholds(self):
        """
        Simulate s residual values and compute s thresholds for each point of a pre-defined grid

        Parameters
        ----------

        Returns
        -------
        """

        # get grid points and evaluated values in the intrusion frame
        grid_points = self.simulation_grid
        grid_points_coord0 = self.simulation_grid[1]
        grid_points_coord1 = self.simulation_grid[2]
        grid_points_coord2 = self.simulation_grid[3]

        # generate data frame containing input data for simulation
        data_sides = self.lateral_contact_data[0]

        minP = min(
            self.vertical_contact_data[0]['coord1'].min(),
            self.vertical_contact_data[1]['coord1'].min(),
            self.lateral_contact_data[0]['coord1'].min())
        maxP = max(
            self.vertical_contact_data[0]['coord1'].max(),
            self.vertical_contact_data[1]['coord1'].max(),
            self.lateral_contact_data[0]['coord1'].max())
        minL = min(
            self.vertical_contact_data[0]['coord2'].min(),
            self.vertical_contact_data[1]['coord2'].min(),
            self.lateral_contact_data[0]['coord2'].min())
        maxL = max(
            self.vertical_contact_data[0]['coord2'].max(),
            self.vertical_contact_data[1]['coord2'].max(),
            self.lateral_contact_data[0]['coord2'].max())

        # -- Min side (s<0)
        data_minL = self.lateral_contact_data[1]
        data_conceptual_minL = self.lateral_extent_model(
            data_minL, minP=minP, maxP=maxP, minS=minL, maxS=maxL
        )
        data_residual_minL = (
            data_conceptual_minL[:, 1] - data_minL.loc[:, "coord2"]
        ).to_numpy()
        inputsimdata_minL = data_minL.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()
        inputsimdata_minL.loc[:, "l_residual"] = data_residual_minL
        inputsimdata_minL.loc[:, "l_conceptual"] = data_conceptual_minL[:, 1]
        inputsimdata_minL.loc[:, "ref_coord"] = 0

        # -- Max side (s>0)
        data_maxL = self.lateral_contact_data[2]
        data_conceptual_maxL = self.lateral_extent_model(
            data_maxL, minP=minP, maxP=maxP, minS=minL, maxS=maxL
        )
        data_residual_maxL = (
            data_conceptual_maxL[:, 0] - data_maxL.loc[:, "coord2"]
        ).to_numpy()
        inputsimdata_maxL = data_maxL.loc[
            :, ["X", "Y", "Z", "coord0", "coord1", "coord2"]
        ].copy()
        inputsimdata_maxL.loc[:, "l_residual"] = data_residual_maxL
        inputsimdata_maxL.loc[:, "l_conceptual"] = data_conceptual_maxL[:, 0]
        inputsimdata_maxL.loc[:, "ref_coord"] = 0

        self.lateral_sgs_input_data = [inputsimdata_minL, inputsimdata_maxL]

        # Compute simulation parameters if not defined
        if self.lateral_sgs_parameters.get("zmin") == None:
            self.lateral_sgs_parameters["zmin"] = min(
                inputsimdata_maxL["l_residual"].min(),
                inputsimdata_minL["l_residual"].min(),
            )
        if self.lateral_sgs_parameters.get("zmax") == None:
            self.lateral_sgs_parameters["zmax"] = max(
                inputsimdata_maxL["l_residual"].max(),
                inputsimdata_minL["l_residual"].max(),
            )
        if self.lateral_sgs_parameters.get("xmn") == None:
            self.lateral_sgs_parameters["xmn"] = np.nanmin(grid_points_coord1)
        if self.lateral_sgs_parameters.get("xsiz") == None:
            nx = self.lateral_sgs_parameters.get("nx")
            minC0 = np.nanmin(grid_points_coord1)
            maxC0 = np.nanmax(grid_points_coord1)
            self.lateral_sgs_parameters["xsiz"] = (maxC0 - minC0) / nx

        # Simulation of lateral extent
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
        zmin = self.lateral_sgs_parameters.get("zmin")
        zmax = self.lateral_sgs_parameters.get("zmax")
        nxdis = self.lateral_sgs_parameters.get("nxdis")
        nydis = self.lateral_sgs_parameters.get("nydis")
        ndmin = self.lateral_sgs_parameters.get("ndmin")
        ndmax = self.lateral_sgs_parameters.get("ndmax")
        radius = self.lateral_sgs_parameters.get("radius")

        if geostats is None:
            raise Exception("geostats is not installed, pip install geostatspy")
        l_min_simulation = geostats.sgsim(
            inputsimdata_minL,
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

        l_max_simulation = geostats.sgsim(
            inputsimdata_maxL,
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

        # self.simulationGSLIB_s_outcome = [s_min_simulation, s_max_simulation]

        # Create dataframe containing S threshold for each grid point

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
            lateral_thresholds, minP=minP, maxP=maxP, minS=minL, maxS=maxL
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
            lateral_thresholds.coord1 < minP,["min_l_threshold", "max_l_threshold"]
            ] = [0.00001, 0.00001]

        lateral_thresholds.loc[
            lateral_thresholds.coord1 > maxP,["min_l_threshold", "max_l_threshold"]
            ] = [0.00001, 0.00001]

        self.lateral_simulated_thresholds = lateral_thresholds

    def simulate_vertical_thresholds(self):
        """
        Simulate g residual values and compute g thresholds for each point of a pre-defined grid.
        Computes two sets of g thresholds: one to constraint the contact opposite to the intrusion network (using residual values),
        and another one to better condition the intrusion network contact to the data.

        Parameters
        ----------

        Returns
        -------
        """

        # Get grid points and evaluated values in the intrusion frame
        grid_points = self.simulation_grid
        grid_points_coord0 = self.simulation_grid[1]
        grid_points_coord1 = self.simulation_grid[2]
        grid_points_coord2 = self.simulation_grid[3]

        # Generate data frame containing input data for simulation
        inet_data = self.vertical_contact_data[0]
        other_contact_data = self.vertical_contact_data[1]

        # --- parameters for conceptual model
        meanG = other_contact_data.loc[:, "coord0"].mean()
        minP = min(inet_data["coord1"].min(), other_contact_data["coord1"].min())
        maxP = max(inet_data["coord1"].max(), other_contact_data["coord1"].max())
        minS = min(inet_data["coord2"].min(), other_contact_data["coord2"].min())
        maxS = max(inet_data["coord2"].max(), other_contact_data["coord2"].max())
        maxG = other_contact_data["coord0"].max()
        coordPS_of_maxG = (
            other_contact_data[
                other_contact_data.coord0 == other_contact_data.coord0.max()
            ]
            .loc[:, ["coord1", "coord2"]]
            .to_numpy()
        )
        vertex = [coordPS_of_maxG[0][0], coordPS_of_maxG[0][1], maxG]

        # --- growth simulation input data (max G, simulation of contact opposite to intrusion network)

        data_conceptual_G = self.vertical_extent_model(
            other_contact_data,
            mean_growth=meanG,
            minP=minP,
            maxP=maxP,
            minS=minS,
            maxS=maxS,
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

        self.simulation_g_inputdata = [inputsimdata_maxG, inputsimdata_inetG]

        # Simulation
        # --- compute simulation parameters if not defined
        if self.simulation_g_parameters.get("nx") == None:
            self.simulation_g_parameters["nx"] = grid_points[4][
                0
            ]  
        if self.simulation_g_parameters.get("ny") == None:
            self.simulation_g_parameters["ny"] = grid_points[4][
                1
            ]  
        if self.simulation_g_parameters.get("xmn") == None:
            self.simulation_g_parameters["xmn"] = np.nanmin(grid_points_coord1)
        if self.simulation_g_parameters.get("ymn") == None:
            self.simulation_g_parameters["ymn"] = np.nanmin(grid_points_coord2)
        if self.simulation_g_parameters.get("xsiz") == None:
            nx = self.simulation_g_parameters.get("nx")
            minC1 = np.nanmin(grid_points_coord1)
            maxC1 = np.nanmax(grid_points_coord1)
            self.simulation_g_parameters["xsiz"] = (maxC1 - minC1) / nx
        if self.simulation_g_parameters.get("ysiz") == None:
            ny = self.simulation_g_parameters.get("ny")
            minC2 = np.nanmin(grid_points_coord2)
            maxC2 = np.nanmax(grid_points_coord2)
            self.simulation_g_parameters["ysiz"] = (maxC2 - minC2) / ny
        if self.simulation_g_parameters.get("zmin") == None:
            self.simulation_g_parameters["zmin"] = inputsimdata_maxG["g_residual"].min()
        if self.simulation_g_parameters.get("zmax") == None:
            self.simulation_g_parameters["zmax"] = inputsimdata_maxG["g_residual"].max()
        if self.simulation_g_parameters.get("zmin2") == None:
            self.simulation_g_parameters["zmin2"] = inputsimdata_inetG["coord0"].min()
        if self.simulation_g_parameters.get("zmax2") == None:
            self.simulation_g_parameters["zmax2"] = inputsimdata_inetG["coord0"].max()

        # --- make variogram
        # self.make_g_simulation_variogram()

        # --- get parameters for simulation
        tmin = self.simulation_g_parameters.get("tmin")
        tmax = self.simulation_g_parameters.get("tmax")
        itrans = self.simulation_g_parameters.get("itrans")
        ktype = self.simulation_g_parameters.get("ktype")
        nx = self.simulation_g_parameters.get("nx")
        ny = self.simulation_g_parameters.get("ny")
        xmn = self.simulation_g_parameters.get("xmn")
        ymn = self.simulation_g_parameters.get("ymn")
        xsiz = self.simulation_g_parameters.get("xsiz")
        ysiz = self.simulation_g_parameters.get("ysiz")
        zmin = self.simulation_g_parameters.get("zmin")
        zmax = self.simulation_g_parameters.get("zmax")
        zmin2 = self.simulation_g_parameters.get("zmin2")
        zmax2 = self.simulation_g_parameters.get("zmax2")
        nxdis = self.simulation_g_parameters.get("nxdis")
        nydis = self.simulation_g_parameters.get("nydis")
        ndmin = self.simulation_g_parameters.get("ndmin")
        ndmax = self.simulation_g_parameters.get("ndmax")
        radius = self.simulation_g_parameters.get("radius")

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
            vario=self.simulation_g_variogram,
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
            vario=self.simulation_g_variogram,
        )

        self.simulationGSLIB_g_outcome = [g_min_simulation, g_max_simulation]
        # Create dataframe containing S threshold for each grid point
        lower_extent_gps = [0, xmn, ymn]
        upper_extent_gps = [0, xmn + xsiz * nx, ymn + ysiz * ny]

        # -- transform SGS output array to Intrusion Frame coordinates
        # -- grid_from_array returns a matrix of [node_i,node_j,g,p,s,values in array]
        g_residual = grid_from_array(
            g_max_simulation, ["X", 0], lower_extent_gps, upper_extent_gps
        )
        g_minimum = grid_from_array(
            g_min_simulation, ["X", 0], lower_extent_gps, upper_extent_gps
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
            minS=minS,
            maxS=maxS,
            vertex=vertex,
        )
        g_maximum = conceptual_maxg[:, 1] - g_residual[:, 5]
        simulation_g_threshold.loc[:, "g_maximum"] = g_maximum

        self.simulated_g_thresholds = simulation_g_threshold


    def evaluate_value(self, points):

        """
        points : array (x,y,z),  points where the intrusion is evaluated.

        """

        simulation_g_data = self.simulation_growth_data
        simulation_s_data = self.lateral_simulated_thresholds
        # inet_polygon = None
        # intrusion_frame = self.intrusion_frame
        model = self.model

        # ---> returns indicator function and scalar field with isovalue 0 = intrusion boundary

        # compute coordinates values for each evaluated point
        intrusion_coord0_pts = self.intrusion_frame[0].evaluate_value(points)
        intrusion_coord1_pts = self.intrusion_frame[1].evaluate_value(points)
        intrusion_coord2_pts = self.intrusion_frame[2].evaluate_value(points)

        # ------ lateral extent thresholds for each of the evaluated points -------------

        # if lateral extent values were simulated
        if simulation_s_data is None:
            print("No simultion for lateral extent")
        else:
            print("Assigning lateral thresholds")
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
            print("No simulation for vertical extent")
        else:
            # containers for thresholds
            print("Assigning vertical thresholds")
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

        mid_point = g_minside_threshold + ((g_maxside_threshold - g_minside_threshold)/2)

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
