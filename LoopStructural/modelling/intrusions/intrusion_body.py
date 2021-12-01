import numpy as np
import pandas as pd

from LoopStructural.utils import getLogger
logger = getLogger(__name__)
#GSLIB library        
try:
    import geostatspy.GSLIB as GSLIB                          # GSLIB utilities, viz and wrapped functions
    import geostatspy.geostats as geostats                    # GSLIB converted to Python

except ImportError:
    logger.error("GeostatPy not installed \n"
                "pip install geostatspy")


#Geometrical conceptual models functions
from LoopStructural.modelling.intrusions.intrusion_support_functions import *

class IntrusionBody():
    def __init__(self, feature_data, name,
                 intrusion_network = None,
                 intrusion_frame = None,
                 model = None):
        """
        Simulates the intrusion lateral and vertical extent, 
        the intrusion structural frame and conceptual geometrical models.

        Parameters
        ----------
        feature_data : dataframe, contains data point that are used to constrain conceptual models
        intrusion_network : IntrusionNetwork object
            intrusion network and its related data
        intrusion_frame : IntrusionFrame
            the intrusion frame defining the coarse-scale intrusion geometry
        """
        
        self.data = feature_data
        self.name = name
        self.intrusion_network = intrusion_network
        self.intrusion_frame = intrusion_frame
        self.model = model
        # conceptual models
        self.lateral_extent_model = None
        self.vertical_extent_model = None
        # containers for data related to each simulation
        self.lateral_contact_data = []
        self.vertical_contact_data = []
        #grid for simulation
        self.simulation_grid = None
        # parameters for each simulation and container for simulated values within the model
        self.simulation_s_parameters = None #input parameters for simulation
        self.simulation_s_variogram = None #inout variogram
        self.simulation_sdata = None  # list [data of all sides, data of side S<0, data of side S>0]
        self.simulation_s_inputdata = [] #residual values
        self.simulationGSLIB_s_outcome = []
        self.simulated_s_thresholds = None #compute thresholds
        
        self.simulation_g_parameters = None
        self.simulation_g_variogram = None
        self.simulation_gdata = None # list [all data, data of intrusion network contact, data of other contact]
        self.simulation_g_inputdata = []
        self.simulationGSLIB_g_outcome = []
        self.simulated_g_thresholds = None
        
    def set_lateral_extent_conceptual_model(self, function = None):
        """Set geometrical conceptual model function to simulate lateral extent
        Function must be an input for the simulation

        Parameters
        ----------
        function : function name
        
        Returns
        ----------
        """
        if function == None:
            logger.error("Specify function of lateral extent conceptual model")
        else:   
            self.lateral_extent_model = function
            
    def set_vertical_extent_conceptual_model(self, function = None):
        """Set geometrical conceptual model function to simulate vertical extent
        Function must be an input for the simulation

        Parameters
        ----------
        function : function name
        
        Returns
        ----------
        """
        if function == None:
            logger.error("Specify function of vertical extent conceptual model")
        else:   
            self.vertical_extent_model = function
    
    def set_data_for_s_simulation(self):
        """Set data for lateral extent (distances in L axis) simulation
        Evaluate data points in the intrusion frame coodinates, 
        and separates data between sides (> or <0)

        Parameters
        ----------
        
        Returns
        ----------
        """
        import datetime

        data = self.data.copy()
        data_xyz = data.loc[:,['X','Y','Z']].to_numpy()
        # data['coord0'] = self.intrusion_frame[0].evaluate_value(self.model.scale(data_xyz, inplace = False))
        # data['coord1'] = self.intrusion_frame[1].evaluate_value(self.model.scale(data_xyz, inplace = False))
        # data['coord2'] = self.intrusion_frame[2].evaluate_value(self.model.scale(data_xyz, inplace = False))
        data.loc[:,'coord0'] = self.intrusion_frame[0].evaluate_value(data_xyz)
        data.loc[:,'coord1'] = self.intrusion_frame[1].evaluate_value(data_xyz)
        data.loc[:,'coord2'] = self.intrusion_frame[2].evaluate_value(data_xyz)
        data_minside = data[(data['intrusion_side'] == 'yes') & (data['coord2']<=0)].copy()
        data_minside.reset_index(inplace = True)

        data_maxside = data[(data['intrusion_side'] == 'yes') & (data['coord2']>0)].copy()
        data_maxside.reset_index(inplace = True)
        
        data_sides = data_minside.append(data_maxside)
        data_sides.reset_index(inplace = True)
        
        self.lateral_contact_data = [data_sides, data_minside, data_maxside]

        return data_sides, data_minside, data_maxside
    
    def set_data_for_g_simulation(self):
        """Set data for vertical extent (distances in G axis) simulation
        Evaluate data points in the intrusion frame coodinates,
        and separate data between intrusion network contact data and other contact data.

        Parameters
        ----------
        
        Returns
        ----------
        """
        intrusion_network_data_xyz = self.intrusion_network.intrusion_network_data.loc[:,['X','Y','Z']].to_numpy()
        intrusion_network_data = self.intrusion_network.intrusion_network_data.loc[:,['X','Y','Z']].copy()
        intrusion_network_data.loc[:,'coord0'] = self.intrusion_frame[0].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:,'coord1'] = self.intrusion_frame[1].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.loc[:,'coord2'] = self.intrusion_frame[2].evaluate_value(intrusion_network_data_xyz)
        intrusion_network_data.reset_index(inplace=True)
        
        other_contact_data_xyz = self.intrusion_network.other_contact_data.loc[:,['X','Y','Z']].to_numpy()
        other_contact_data = self.intrusion_network.other_contact_data.loc[:,['X','Y','Z']].copy()
        other_contact_data.loc[:,'coord0'] = self.intrusion_frame[0].evaluate_value(other_contact_data_xyz)
        other_contact_data.loc[:,'coord1'] = self.intrusion_frame[1].evaluate_value(other_contact_data_xyz)
        other_contact_data.loc[:,'coord2'] = self.intrusion_frame[2].evaluate_value(other_contact_data_xyz)
        other_contact_data.reset_index(inplace=True)
        
        self.vertical_contact_data = [intrusion_network_data, other_contact_data]  

            
    def create_grid_for_simulation(self, spacing = None):
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
            
        x = np.linspace(self.model.origin[0], self.model.maximum[0], spacing[0])
        y = np.linspace(self.model.origin[1], self.model.maximum[1], spacing[1])
        z = np.linspace(self.model.origin[2], self.model.maximum[2], spacing[2])

        xx,yy,zz = np.meshgrid(x,y,z)
        grid_points = np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T
        
        grid_points_coord0 = self.intrusion_frame[0].evaluate_value(self.model.scale(grid_points, inplace = False))
        grid_points_coord1 = self.intrusion_frame[1].evaluate_value(self.model.scale(grid_points, inplace = False))
        grid_points_coord2 = self.intrusion_frame[2].evaluate_value(self.model.scale(grid_points, inplace = False))
        
        self.simulation_grid = [grid_points, grid_points_coord0, grid_points_coord1, grid_points_coord2, spacing]

    
    def set_s_simulation_GSLIBparameters(self, lateral_simulation_parameters):        
       
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
        if 'tmin' in lateral_simulation_parameters: tmin = lateral_simulation_parameters['tmin']
        else: tmin = -999
        
        if 'tmax' in lateral_simulation_parameters: tmax = lateral_simulation_parameters['tmax']
        else: tmax = 999

        if 'itrans' in lateral_simulation_parameters: itrans = lateral_simulation_parameters['itrans']
        else: itrans = 1

        if 'ktype' in lateral_simulation_parameters: ktype = lateral_simulation_parameters['ktype']
        else: ktype = 0
        
        if 'nx' in lateral_simulation_parameters: nx = lateral_simulation_parameters['nx']
        else: nx = 200
        
        if 'ny' in lateral_simulation_parameters: ny= lateral_simulation_parameters['ny']
        else: ny = 3

        if 'xmn' in lateral_simulation_parameters: xmn = lateral_simulation_parameters['xmn']
        else: xmn = None

        if 'ymn' in lateral_simulation_parameters: ymn = lateral_simulation_parameters['ymn']
        else: ymn = -0.0001

        if 'xsiz' in lateral_simulation_parameters: xsiz = lateral_simulation_parameters['xsiz']
        else: xsiz = None
        
        if 'ysiz' in lateral_simulation_parameters: ysiz = lateral_simulation_parameters['ysiz']
        else: ysiz = 0.0001

        if 'zmin' in lateral_simulation_parameters: zmin = lateral_simulation_parameters['zmin']
        else: zmin = None

        if 'zmax' in lateral_simulation_parameters: zmax = lateral_simulation_parameters['zmax']
        else: zmax = None

        if 'nxdis' in lateral_simulation_parameters: nxdis = lateral_simulation_parameters['nxdis']
        else: nxdis = 1

        if 'nydis' in lateral_simulation_parameters: nydis = lateral_simulation_parameters['nydis']
        else: nydis = 1

        if 'ndmin' in lateral_simulation_parameters: ndmin = lateral_simulation_parameters['ndmin']
        else: ndmin = 0

        if 'ndmax' in lateral_simulation_parameters: ndmax = lateral_simulation_parameters['ndmax']
        else: ndmax = 3

        if 'radius' in lateral_simulation_parameters: radius= lateral_simulation_parameters['radius']
        else: radius= 500


        s_parameters = {'tmin':tmin , 'tmax':tmax, 'itrans':itrans, 'ktype':ktype, 
                        'nx':nx, 'ny':ny, 'xmn':xmn, 'ymn':ymn, 'xsiz' :xsiz, 'ysiz': ysiz,
                        'zmin':zmin, 'zmax':zmax, 'nxdis':nxdis, 'nydis':nydis,
                        'ndmin':ndmin, 'ndmax':ndmax, 'radius':radius}
        
        self.simulation_s_parameters = s_parameters
        
        
    def set_g_simulation_GSLIBparameters(self, vertical_simulation_parameters):
    
        
        """
        Simulation parameters for vertical extent simulation

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
        if 'tmin' in vertical_simulation_parameters: tmin = vertical_simulation_parameters['tmin']
        else: tmin = -999
        
        if 'tmax' in vertical_simulation_parameters: tmax = vertical_simulation_parameters['tmax']
        else: tmax = 999

        if 'itrans' in vertical_simulation_parameters: itrans = vertical_simulation_parameters['itrans']
        else: itrans = 1

        if 'ktype' in vertical_simulation_parameters: ktype = vertical_simulation_parameters['ktype']
        else: ktype = 0
        
        if 'nx' in vertical_simulation_parameters: nx = vertical_simulation_parameters['nx']
        else: nx = None
        
        if 'ny' in vertical_simulation_parameters: ny= vertical_simulation_parameters['ny']
        else: ny = None

        if 'xmn' in vertical_simulation_parameters: xmn = vertical_simulation_parameters['xmn']
        else: xmn = None

        if 'ymn' in vertical_simulation_parameters: ymn = vertical_simulation_parameters['ymn']
        else: ymn = None

        if 'xsiz' in vertical_simulation_parameters: xsiz = vertical_simulation_parameters['xsiz']
        else: xsiz = None
        
        if 'ysiz' in vertical_simulation_parameters: ysiz = vertical_simulation_parameters['ysiz']
        else: ysiz = None

        if 'zmin' in vertical_simulation_parameters: zmin = vertical_simulation_parameters['zmin']
        else: zmin = None

        if 'zmax' in vertical_simulation_parameters: zmax = vertical_simulation_parameters['zmax']
        else: zmax = None

        if 'zmin2' in vertical_simulation_parameters: zmin2 = vertical_simulation_parameters['zmin2']
        else: zmin2 = None

        if 'zmax2' in vertical_simulation_parameters: zmax2 = vertical_simulation_parameters['zmax2']
        else: zmax2 = None

        if 'nxdis' in vertical_simulation_parameters: xdis = vertical_simulation_parameters['nxdis']
        else: nxdis = 1

        if 'nydis' in vertical_simulation_parameters: nydis = vertical_simulation_parameters['nydis']
        else: nydis = 1

        if 'ndmin' in vertical_simulation_parameters: ndmin = vertical_simulation_parameters['ndmin']
        else: ndmin = 0

        if 'ndmax' in vertical_simulation_parameters: ndmax = vertical_simulation_parameters['ndmax']
        else: ndmax = 3

        if 'radius' in vertical_simulation_parameters: radius= vertical_simulation_parameters['radius']
        else: radius= 500
        
        g_parameters = {'tmin':tmin , 'tmax':tmax, 'itrans':itrans, 'ktype':ktype, 
                        'nx':nx, 'ny':ny, 'xmn':xmn, 'ymn':ymn, 'xsiz' :xsiz, 'ysiz': ysiz,
                        'zmin':zmin, 'zmax':zmax, 'zmin2':zmin2, 'zmax2':zmax2, 'nxdis':nxdis, 'nydis':nydis,
                        'ndmin':ndmin, 'ndmax':ndmax, 'radius':radius}
        
        self.simulation_g_parameters = g_parameters
    
    def make_s_simulation_variogram(self, lateral_simulation_parameters):
        """
        Make variogram for lateral extent simulation
        By default: variogram with no nugget effect, 1 spherical nested structure, isotropic and infinite range

        Parameters
        ----------
        python dictionary with the following parameters -->

        nug = 0; nst = 1                   # nugget effect = 0 and 1 nested structure
        it1 = 1                            # nested structure type (1 - spherical, 2 - exponential, 3 - Gaussian)
        cc1 = 1                            # One structure and no nugget 
        azi1 = 90                          # azimuth of this nested structure
        hmaj1 = 9999 ; hmin1 = 9999   # range of this nested structure in the major (hmaj1) and minor (hmin1) direction
        
        Returns
        ----
        """ 
        if 'nug' in lateral_simulation_parameters: nug = lateral_simulation_parameters['nug']
        else: nug = 0

        if 'nst' in lateral_simulation_parameters: nst = lateral_simulation_parameters['nst']
        else: nst = 1

        if 'it1' in lateral_simulation_parameters: it1 = lateral_simulation_parameters['it1']
        else: it1 = 1

        if 'cc1' in lateral_simulation_parameters: cc1 = lateral_simulation_parameters['cc1']
        else: cc1 = 1

        if 'azi1' in lateral_simulation_parameters: azi1 = lateral_simulation_parameters['azi1']
        else: azi1 = 90

        if 'hmaj1' in lateral_simulation_parameters: hmaj1 = lateral_simulation_parameters['hmaj1']
        else: hmaj1 = 999999

        if 'hmin1' in lateral_simulation_parameters: hmin1 = lateral_simulation_parameters['hmin1']
        else: hmin1 = 999999
            

        variogram = GSLIB.make_variogram(nug,nst,it1,cc1,azi1,hmaj1,hmin1)
        self.simulation_s_variogram = variogram
        
    def make_g_simulation_variogram(self, vertical_simulation_parameters):
        """
        Make variogram for vertical extent simulation
        By default: variogram with no nugget effect, 1 spherical nested structure, isotropic and infinite range

        Parameters
        ----------
        nug = 0; nst = 1                   # nugget effect = 0 and 1 nested structure
        it1 = 1                            # nested structure type (1 - spherical, 2 - exponential, 3 - Gaussian)
        cc1 = 1                            # One structure and no nugget 
        azi1 = 90                          # azimuth of this nested structure
        hmaj1 = 9999 ; hmin1 = 9999   # range of this nested structure in the major (hmaj1) and minor (hmin1) direction
        
        Returns
        ----
        """  

        if 'nug' in vertical_simulation_parameters: nug = vertical_simulation_parameters['nug']
        else: nug = 0

        if 'nst' in vertical_simulation_parameters: nst = vertical_simulation_parameters['nst']
        else: nst = 1

        if 'it1' in vertical_simulation_parameters: it1 = vertical_simulation_parameters['it1']
        else: it1 = 1

        if 'cc1' in vertical_simulation_parameters: cc1 = vertical_simulation_parameters['cc1']
        else: cc1 = 1

        if 'azi1' in vertical_simulation_parameters: azi1 = vertical_simulation_parameters['azi1']
        else: azi1 = 90

        if 'hmaj1' in vertical_simulation_parameters: hmaj1 = vertical_simulation_parameters['hmaj1']
        else: hmaj1 = 999999

        if 'hmin1' in vertical_simulation_parameters: hmin1 = vertical_simulation_parameters['hmin1']
        else: hmin1 = 999999    
        
        variogram = GSLIB.make_variogram(nug,nst,it1,cc1,azi1,hmaj1,hmin1)
        self.simulation_g_variogram = variogram
    
        
    def simulate_s_thresholds(self):
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
        minP = data_sides['coord1'].min() ;  maxP = data_sides['coord1'].max()
        minS = data_sides['coord2'].min() ;  maxS = data_sides['coord2'].max()
        
        # -- Min side (s<0)
        data_minS = self.lateral_contact_data[1]
        data_conceptual_minS = self.lateral_extent_model(data_minS, minP=minP, maxP=maxP, minS=minS, maxS=maxS)
        data_residual_minS = (data_conceptual_minS[:,1] - data_minS.loc[:,'coord2']).to_numpy()
        inputsimdata_minS = data_minS.loc[:, ['X','Y','Z','coord0','coord1','coord2']].copy()
        inputsimdata_minS.loc[:,'s_residual'] = data_residual_minS
        inputsimdata_minS.loc[:,'s_conceptual'] = data_conceptual_minS[:,1]
        inputsimdata_minS.loc[:,'ref_coord'] = 0               
        
        # -- Max side (s>0)
        data_maxS = self.lateral_contact_data[2]
        data_conceptual_maxS = self.lateral_extent_model(data_maxS, minP=minP, maxP=maxP, minS=minS, maxS=maxS)
        data_residual_maxS = (data_conceptual_maxS[:,0] - data_maxS.loc[:,'coord2']).to_numpy()
        inputsimdata_maxS = data_maxS.loc[:, ['X','Y','Z','coord0','coord1','coord2']].copy()
        inputsimdata_maxS.loc[:,'s_residual'] = data_residual_maxS
        inputsimdata_maxS.loc[:,'s_conceptual'] = data_conceptual_maxS[:,0]
        inputsimdata_maxS.loc[:,'ref_coord'] = 0
        
        self.simulation_s_inputdata = [inputsimdata_minS, inputsimdata_maxS]
        
        # Compute simulation parameters if not defined
        if self.simulation_s_parameters.get('zmin') == None:
            self.simulation_s_parameters['zmin'] = min(inputsimdata_maxS['s_residual'].min(),inputsimdata_minS['s_residual'].min())
        if self.simulation_s_parameters.get('zmax') == None:
            self.simulation_s_parameters['zmax']= max(inputsimdata_maxS['s_residual'].max(),inputsimdata_minS['s_residual'].max())
        if self.simulation_s_parameters.get('xmn') == None:
            self.simulation_s_parameters['xmn'] = np.nanmin(grid_points_coord1)
        if self.simulation_s_parameters.get('xsiz') == None:
            nx = self.simulation_s_parameters.get('nx')
            minC0 = np.nanmin(grid_points_coord1)
            maxC0 = np.nanmax(grid_points_coord1)
            self.simulation_s_parameters['xsiz'] = (maxC0-minC0)/nx
            
        # Make variogram
        # self.make_s_simulation_variogram()
        
        # Simulation of lateral extent
        tmin = self.simulation_s_parameters.get('tmin'); tmax = self.simulation_s_parameters.get('tmax')
        itrans = self.simulation_s_parameters.get('itrans'); ktype = self.simulation_s_parameters.get('ktype')                    
        nx = self.simulation_s_parameters.get('nx') ; ny = self.simulation_s_parameters.get('ny')
        xmn = self.simulation_s_parameters.get('xmn') ; ymn = self.simulation_s_parameters.get('ymn')
        xsiz = self.simulation_s_parameters.get('xsiz') ; ysiz = self.simulation_s_parameters.get('ysiz')
        zmin = self.simulation_s_parameters.get('zmin'); zmax = self.simulation_s_parameters.get('zmax')
        nxdis = self.simulation_s_parameters.get('nxdis'); nydis = self.simulation_s_parameters.get('nydis')
        ndmin = self.simulation_s_parameters.get('ndmin'); ndmax = self.simulation_s_parameters.get('ndmax')         
        radius = self.simulation_s_parameters.get('radius') 
        
        if geostats is None:
            raise Exception('geostats is not installed')    
        s_min_simulation = geostats.sgsim(inputsimdata_minS,'coord1','ref_coord','s_residual',wcol=-1,scol=-1,
                                          tmin=tmin,tmax=tmax,itrans=itrans,ismooth=0,dftrans=0,tcol=0,twtcol=0,
                                          zmin=zmin,zmax=zmax,ltail=1,ltpar=0.0,utail=1,utpar=0.3,nsim=1,
                                          nx=nx,xmn=xmn,xsiz=xsiz,ny=ny,ymn=ymn,ysiz=ysiz,
                                          seed=73073,
                                          ndmin=ndmin,ndmax=ndmax,nodmax=1,mults=0,nmult=2,noct=-1,radius=radius,radius1=10,sang1=0,
                                          mxctx=1,mxcty=1,ktype=ktype,colocorr=0.0,sec_map=0,vario=self.simulation_s_variogram)
        
        s_max_simulation = geostats.sgsim(inputsimdata_maxS,'coord1','ref_coord','s_residual',wcol=-1,scol=-1,
                                          tmin=tmin,tmax=tmax,itrans=itrans,ismooth=0,dftrans=0,tcol=0,twtcol=0,
                                          zmin=zmin,zmax=zmax,ltail=1,ltpar=0.0,utail=1,utpar=0.3,nsim=1,
                                          nx=nx,xmn=xmn,xsiz=xsiz,ny=ny,ymn=ymn,ysiz=ysiz,
                                          seed=73073,
                                          ndmin=ndmin,ndmax=ndmax,nodmax=1,mults=0,nmult=2,noct=-1,radius=radius,radius1=10,sang1=0,
                                          mxctx=1,mxcty=1,ktype=ktype,colocorr=0.0,sec_map=0,vario=self.simulation_s_variogram)

        self.simulationGSLIB_s_outcome = [s_min_simulation, s_max_simulation]
        
        # Create dataframe containing S threshold for each grid point

        propagation_grid_model = np.linspace(xmn,xmn + (nx*xsiz), nx)

        simulation_s_thresholds = pd.DataFrame(columns = ['coord1','min_s_residual','min_s_threshold','max_s_residual','max_s_threshold'])
        simulation_s_thresholds['coord1'] = propagation_grid_model

        model_conceptual_s = self.lateral_extent_model(simulation_s_thresholds, minP=minP, maxP=maxP, minS=minS, maxS=maxS)
        simulation_s_thresholds['conceptual_mins'] = model_conceptual_s[:,1]
        simulation_s_thresholds['min_s_residual'] = s_min_simulation[1]
        simulation_s_thresholds['min_s_threshold'] = model_conceptual_s[:,1] - s_min_simulation[1]
        
        simulation_s_thresholds['conceptual_maxs'] = model_conceptual_s[:,0]
        simulation_s_thresholds['max_s_residual'] = s_max_simulation[1]
        simulation_s_thresholds['max_s_threshold'] = model_conceptual_s[:,0] - s_max_simulation[1]
        simulation_s_thresholds.sort_values(['coord1'], ascending = [True], inplace = True)
        
        # Ignore simulated data outside area covered by input data
        for j in range(len(simulation_s_thresholds)):
            if simulation_s_thresholds.loc[j,'coord1'] < minP:
                simulation_s_thresholds.loc[j, ['min_s_threshold','max_s_threshold']] = [0.00001,0.00001]   
            if simulation_s_thresholds.loc[j,'coord1'] > maxP:
                simulation_s_thresholds.loc[j, ['min_s_threshold','max_s_threshold']] = [0.00001,0.00001]
        
        self.simulated_s_thresholds = simulation_s_thresholds
        
#         return simulation_s_thresholds
    
    def simulate_g_thresholds(self):
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
        meanG = other_contact_data.loc[:,'coord0'].mean()
        minP = min(inet_data['coord1'].min(), other_contact_data['coord1'].min())
        maxP = max(inet_data['coord1'].max(), other_contact_data['coord1'].max())
        minS = min(inet_data['coord2'].min(), other_contact_data['coord2'].min())
        maxS = max(inet_data['coord2'].max(), other_contact_data['coord2'].max())
        maxG =  other_contact_data['coord0'].max()
        coordPS_of_maxG = other_contact_data[other_contact_data.coord0 == other_contact_data.coord0.max()].loc[:,['coord1','coord2']].to_numpy()
        vertex = [coordPS_of_maxG[0][0],coordPS_of_maxG[0][1],maxG]
        
        # --- growth simulation input data (max G, simulation of contact opposite to intrusion network)
        # print(other_contact_data)
        # print(minP,maxP,minS,maxS,vertex)

        data_conceptual_G = self.vertical_extent_model(other_contact_data, mean_growth=meanG, minP=minP, maxP=maxP, minS=minS, maxS=maxS, vertex=vertex)
        data_residual_G = (data_conceptual_G[:,1] - other_contact_data.loc[:,'coord0']).to_numpy()
        inputsimdata_maxG = other_contact_data.loc[:, ['X','Y','Z','coord0','coord1','coord2']].copy()
        inputsimdata_maxG.loc[:,'g_residual'] = data_residual_G
        inputsimdata_maxG.loc[:,'g_conceptual'] = data_conceptual_G[:,1]         
        
        # --- growth simulation input data for intrusion network conditioning
        inputsimdata_inetG = inet_data.loc[:, ['X','Y','Z','coord0','coord1','coord2']].copy()
        
        self.simulation_g_inputdata = [inputsimdata_maxG, inputsimdata_inetG]
    
        # Simulation
        # --- compute simulation parameters if not defined
        if self.simulation_g_parameters.get('nx') == None:
            self.simulation_g_parameters['nx'] = grid_points[4][0] #*2 # grid X spacing
        if self.simulation_g_parameters.get('ny') == None:
            self.simulation_g_parameters['ny'] = grid_points[4][1] #*2 # grid Y spacing
        if self.simulation_g_parameters.get('xmn') == None:
            self.simulation_g_parameters['xmn'] = np.nanmin(grid_points_coord1)
        if self.simulation_g_parameters.get('ymn') == None:
            self.simulation_g_parameters['ymn'] = np.nanmin(grid_points_coord2)
        if self.simulation_g_parameters.get('xsiz') == None:
            nx = self.simulation_g_parameters.get('nx')
            minC1 = np.nanmin(grid_points_coord1)
            maxC1 = np.nanmax(grid_points_coord1)
            self.simulation_g_parameters['xsiz'] = (maxC1-minC1)/nx         
        if self.simulation_g_parameters.get('ysiz') == None:
            yx = self.simulation_g_parameters.get('yx')
            minC2 = np.nanmin(grid_points_coord2)
            maxC2 = np.nanmax(grid_points_coord2)
            self.simulation_g_parameters['ysiz'] = (maxC2-minC2)/nx    
        if self.simulation_g_parameters.get('zmin') == None:
            self.simulation_g_parameters['zmin'] = inputsimdata_maxG['g_residual'].min()
        if self.simulation_g_parameters.get('zmax') == None:
            self.simulation_g_parameters['zmax']= inputsimdata_maxG['g_residual'].max()
        if self.simulation_g_parameters.get('zmin2') == None:
            self.simulation_g_parameters['zmin2'] = inputsimdata_inetG['coord0'].min()
        if self.simulation_g_parameters.get('zmax2') == None:
            self.simulation_g_parameters['zmax2']= inputsimdata_inetG['coord0'].max()

            
        # --- make variogram
        # self.make_g_simulation_variogram()
        
        # --- get parameters for simulation
        tmin = self.simulation_g_parameters.get('tmin'); tmax = self.simulation_g_parameters.get('tmax')
        itrans = self.simulation_g_parameters.get('itrans'); ktype = self.simulation_g_parameters.get('ktype')                    
        nx = self.simulation_g_parameters.get('nx') ; ny = self.simulation_g_parameters.get('ny')
        xmn = self.simulation_g_parameters.get('xmn') ; ymn = self.simulation_g_parameters.get('ymn')
        xsiz = self.simulation_g_parameters.get('xsiz') ; ysiz = self.simulation_g_parameters.get('ysiz')
        zmin = self.simulation_g_parameters.get('zmin'); zmax = self.simulation_g_parameters.get('zmax')
        zmin2 = self.simulation_g_parameters.get('zmin2'); zmax2 = self.simulation_g_parameters.get('zmax2')
        nxdis = self.simulation_g_parameters.get('nxdis'); nydis = self.simulation_g_parameters.get('nydis')
        ndmin = self.simulation_g_parameters.get('ndmin'); ndmax = self.simulation_g_parameters.get('ndmax')         
        radius = self.simulation_g_parameters.get('radius') 
        
        # --- simulation of residual values related to other contact( opposite to intrusion network)
        g_max_simulation = geostats.sgsim(inputsimdata_maxG,'coord1','coord2','g_residual',wcol=-1,scol=-1,
                                          tmin=tmin,tmax=tmax,itrans=itrans,ismooth=0,dftrans=0,tcol=0,twtcol=0,
                                          zmin=zmin,zmax=zmax,ltail=1,ltpar=0.0,utail=1,utpar=0.3,nsim=1,
                                          nx=nx,xmn=xmn,xsiz=xsiz,ny=ny,ymn=ymn,ysiz=ysiz,
                                          seed=73073,
                                          ndmin=ndmin,ndmax=ndmax,nodmax=1,mults=0,nmult=2,noct=-1,radius=radius,radius1=10,sang1=0,
                                          mxctx=1,mxcty=1,ktype=ktype,colocorr=0.0,sec_map=0,vario=self.simulation_g_variogram)
        
        # --- simulation to improve conditioning of intrusion network contact
        g_min_simulation = geostats.sgsim(inputsimdata_inetG,'coord1','coord2','coord0',wcol=-1,scol=-1,
                                          tmin=tmin,tmax=tmax,itrans=itrans,ismooth=0,dftrans=0,tcol=0,twtcol=0,
                                          zmin=zmin2,zmax=zmax2,ltail=1,ltpar=0.0,utail=1,utpar=0.3,nsim=1,
                                          nx=nx,xmn=xmn,xsiz=xsiz,ny=ny,ymn=ymn,ysiz=ysiz,
                                          seed=73073,
                                          ndmin=ndmin,ndmax=ndmax,nodmax=1,mults=0,nmult=2,noct=-1,radius=radius,radius1=10,sang1=0,
                                          mxctx=1,mxcty=1,ktype=ktype,colocorr=0.0,sec_map=0,vario=self.simulation_g_variogram)

        self.simulationGSLIB_g_outcome = [g_min_simulation, g_max_simulation]                 
        # Create dataframe containing S threshold for each grid point
        lower_extent_gps = [0, xmn, ymn]
        upper_extent_gps = [0, xmn+xsiz*nx, ymn+ysiz*ny]
        
        # -- transform SGS output array to Intrusion Frame coordinates
        # -- grid_from_array returns a matrix of [node_i,node_j,g,p,s,values in array]
        g_residual = grid_from_array(g_max_simulation, ['X',0], lower_extent_gps, upper_extent_gps)
        g_minimum = grid_from_array(g_min_simulation, ['X',0], lower_extent_gps, upper_extent_gps)
        
        simulation_g_threshold = pd.DataFrame(columns = ['coord0', 'coord1', 'coord2', 'g_residual', 'g_maximum', 'g_minimum'])
        simulation_g_threshold.loc[:,'coord0'] = g_residual[:,2]
        simulation_g_threshold.loc[:,'coord1'] = g_residual[:,3]
        simulation_g_threshold.loc[:,'coord2'] = g_residual[:,4]
        simulation_g_threshold.loc[:,'g_residual'] = g_residual[:,5]
        simulation_g_threshold.loc[:,'g_minimum'] = g_minimum[:,5]
        conceptual_maxg = self.vertical_extent_model(simulation_g_threshold, mean_growth=meanG, minP=minP, maxP=maxP, minS=minS, maxS=maxS, vertex=vertex)
        g_maximum = conceptual_maxg[:,1] - g_residual[:,5]
        simulation_g_threshold.loc[:,'g_maximum'] = g_maximum
        
        self.simulated_g_thresholds = simulation_g_threshold
        
#         return simulation_g_threshold
        

    