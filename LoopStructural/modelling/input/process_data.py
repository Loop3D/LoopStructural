import pandas as pd
import numpy as np
import os
import logging
import networkx
from scipy.stats import truncnorm
import copy
from .fault_network import FaultNetwork
from LoopStructural.utils import strike_dip_vector
from LoopStructural.utils import getLogger
logger = getLogger(__name__)
class ProcessInputData:
    def __init__(   self, 
                    contacts=None, 
                    contact_orientations = None, 
                    stratigraphic_order = None,
                    fault_orientations = None, 
                    fault_locations = None, 
                    fault_properties = None, 
                    fault_edges = None,
                    intrusions = None,
                    fault_stratigraphy = None,
                    thicknesses = None,
                    colours = None,
                    use_thickness = None,
                    fault_edge_properties=None,
                ):
        """Object to generate loopstructural input dataset from a geological map

        Parameters
        ----------
        contacts : DataFrame
            x,y,z,name for each contact
        contact_orientations : DataFrame
            x,y,z,strike,dip,name for each contact  
        stratigraphic_order : nested list
            a nested list e.g. [['a','b','c'],['d','e']]
            a->b->c are the youngest supergroup, d->e older.
        fault_orientations : DataFrame, optional
            data frame with x,y,z,strike,fault_name, by default None
        fault_locations : DataFrame, optional
            data frame with x,y,z,fault_name, by default None
        fault_properties : DataFrame, optional
            dataframe with properties used for building fault
        fault_edges : list, optional
            intersections between faults, by default None
        intrusions : list, optional
            list containing any intrusions, assumes data in contacts/orientations, by default None
        fault_stratigraphy : graph, optional
            fault stratigraphy intersection, by default None
        thicknesses : dict, optional
            dictionary with stratigraphic thicknesses, by default None
        colours : dict, optional
            dictionary with stratigraphic unit colours, by default None
        use_thickness : bool, optional
            assign value to stratigraphic unit or use interface constraints, by default None
        Notes
        ------
        
        The processor will generate the best possible data set given the input data. If you only want to build a fault 
        network then only fault locations, orientations edges and properties are required
        """
        
        self._stratigraphic_order = stratigraphic_order
        self._thicknesses = thicknesses
        self._use_thickness = use_thickness
        if self.thicknesses is None:
            self._use_thickness = False
        if self._use_thickness is None:
            self._use_thickness = True
        self._vector_scale = 1.
        self._gradient = False
        # all of the properties are processed using setter and getter methods
        # using the property decorator. This means that you can changet he underlying data
        # e.g. the contact dataframe and the model input data will be updated.
        self._contacts = None
        self.contacts = contacts
        self._contact_orientations = None
        self.contact_orientations = contact_orientations
        self._fault_orientations = None
        self.fault_orientations = fault_orientations
        self._fault_locations = None
        self.fault_locations = fault_locations
        # self._fault_dimensions = None
        self._fault_properties = None
        self.fault_properties = fault_properties
        self._fault_network = None
        self.set_fault_network(fault_edges,fault_edge_properties)# = fault_graph
        self._fault_stratigraphy = fault_stratigraphy
        self._intrusions = intrusions
        self._thicknesses = thicknesses
        self._data = None
        self._colours = {}
        self.colours = colours
        # flags
        self._foliation_properties = {} #empty dictionary of foliation parameters
        self.foliation_properties = None

    
    @property
    def colours(self):
        return self._colours

    @colours.setter
    def colours(self,colours):
        if colours is None:
            self._colours = {}
            for s in self.stratigraphic_name:
                self._colours[s] = np.random.random(3)
        else:
            self._colours=colours

    @property
    def stratigraphic_column(self):
        stratigraphic_column = {}
        # add stratigraphy into the column
        unit_id = 0
        val = self._stratigraphic_value()
        for name, sg in self._stratigraphic_order:
            stratigraphic_column[name] = {}
            for g in reversed(sg):
                stratigraphic_column[name][g] = {'max': val[g]+self.thicknesses[g], 'min': val[g] , 'id': unit_id, 'colour':self.colours[g]}
                unit_id += 1
        # add faults into the column 
        stratigraphic_column['faults'] = self.fault_properties.to_dict('index')
        return stratigraphic_column

    @property
    def foliation_properties(self):
        return self._foliation_properties

    @foliation_properties.setter
    def foliation_properties(self,foliation_properties):
        if foliation_properties is None:
            for k in self.stratigraphic_column.keys():
                if k != 'faults':
                    self._foliation_properties[k] = {}
        else:
            self._foliation_properties = foliation_properties

    @property
    def fault_properties(self):
        return self._fault_properties

    @fault_properties.setter
    def fault_properties(self, fault_properties):
        if fault_properties is None:
            return
        
        fault_properties = fault_properties.copy()
        if 'centreEasting' not in fault_properties.columns or \
           'centreNorthing' not in fault_properties.columns or \
           'centreAltitude' not in fault_properties:
            fault_properties['centreEasting'] = np.nan
            fault_properties['centreNorthing'] = np.nan
            fault_properties['centreAltitude'] = np.nan
            for fname in fault_properties.index:
                pts = self.fault_locations.loc[self.fault_locations['feature_name'] == fname,
                                                        ['X','Y','Z']]
                fault_properties.loc[fname,['centreEasting','centreNorthing','centreAltitude']] = np.nanmean(pts,axis=0)
        if 'avgNormalEasting' not in fault_properties.columns or \
           'avgNormalNorthing' not in fault_properties.columns or \
           'avgNormalAltitude' not in fault_properties:
            fault_properties['avgNormalEasting'] = np.nan
            fault_properties['avgNormalNorthing'] = np.nan
            fault_properties['avgNormalAltitude'] = np.nan

            for fname in fault_properties.index:
                pts = self.fault_orientations.loc[self.fault_orientations['feature_name'] == fname,
                                                        ['gx','gy','gz']]
                fault_properties.loc[fname,['avgNormalEasting','avgNormalNorthing','avgNormalAltitude']] = np.nanmean(pts,axis=0)
        if 'avgSlipDirEasting' not in fault_properties.columns or \
           'avgSlipDirNorthing' not in fault_properties.columns or \
           'avgSlipDirAltitude' not in fault_properties:
            fault_properties['avgSlipDirEasting'] = np.nan
            fault_properties['avgSlipDirNorthing'] = np.nan
            fault_properties['avgSlipDirAltitude'] = np.nan
            for fname in fault_properties.index:
                if  'sx' not in self.fault_orientations.columns or \
                    'sy' not in self.fault_orientations.columns or \
                    'sz' not in self.fault_orientations:
                    # if we don't know slip assume its down
                    fault_properties.loc[fname,['avgSlipDirEasting','avgSlipDirNorthing','avgSlipDirAltitude']] = np.array([0.,0.,-1.])
                else:
                    pts = self.fault_orientations.loc[self.fault_orientations['feature_name'] == fname,
                                                            ['sx','sy','sz']]
                    fault_properties.loc[fname,['avgSlipDirEasting','avgSlipDirNorthing','avgSlipDirAltitude']] = np.mean(pts,axis=0)
        self._fault_properties = fault_properties
        self._fault_network = FaultNetwork(list(self.fault_properties.index))
    @property
    def fault_network(self):
        return self._fault_network

    def set_fault_network(self, edges,edge_properties=None):
        if self._fault_network is None:
            self._fault_network = FaultNetwork(list(self.fault_properties.index))
        if edge_properties is None:
            edge_properties = [{} for i in range(len(edges))]
        if edges is None:
            return
        # self._fault_network = FaultNetwork(list(fault_network.nodes))
        for i, e in enumerate(edges):
            self._fault_network.add_connection(e[0],e[1],edge_properties[i])

    def fault_interesections_angle(self,fault1,fault2):
        return np.abs(self.fault_properties.loc[fault1,'dip_dir']-self.fault_properties.loc[fault2,'dip_dir'])

    @property
    def fault_names(self):
        return list(self.fault_properties.index)
        
    @property
    def data(self):
        """This is the main function that does all the work, should be called
        before any of the calculated attributes are accessed
        """
        dataframes = []
        dataframes.append(self.contacts)
        dataframes.append(self.contact_orientations)
        if self.fault_orientations is not None:
            dataframes.append(self.fault_orientations)
        if self.fault_locations is not None:
            dataframes.append(self.fault_locations)
        data = pd.concat(dataframes)
        data.reset_index(inplace=True)
        return data

    @property
    def thicknesses(self):
        return self._thicknesses
    
    @thicknesses.setter
    def thicknesses(self,thicknesses):
        self._thicknesses = thicknesses
        

    @property
    def vector_scale(self):
        return self._vector_scale

    @vector_scale.setter
    def vector_scale(self,vector_scale):
        self._vector_scale = vector_scale
  
    @property
    def stratigraphic_name(self):
        names = []
        for sg in self._stratigraphic_order:
            for g in sg:
                names.append(g)
        return names
    
    
    
    def _stratigraphic_value(self):
        """Calculate the stratigraphic values using stratigraphic order and
        thickness

        Returns
        -------
        dict
            keys are unit name, value is cumulative thickness/implicit function value
        """
        stratigraphic_value = {}
        for name, sg in self._stratigraphic_order:
            value = 0 #reset for each supergroup
            for g in reversed(sg):
                stratigraphic_value[g] = value #+ self._thicknesses[g]
                value+=self._thicknesses[g]
        return stratigraphic_value

    def _update_feature_names(self,dataframe):
            """Function for populating the feature name for stratigraphy using the stratigraphic
            order.

            Parameters
            ----------
            dataframe : DataFrame
                the dataframe to add the new column to
            """
            dataframe['feature_name'] = None
            for name, sg in self._stratigraphic_order:
                for g in sg:
                    dataframe.loc[dataframe['name']==g,'feature_name'] = name
    
    @property
    def contacts(self):
        return self._contacts

    @contacts.setter
    def contacts(self,contacts):
        """Function to convert input contact to loopstructural input

        either uses the thickness values or assigns unique ids given
        the units named in stratigraphic order

        Returns
        -------
        DataFrame
            data frame with x,y,y,val/interface,feature_name
        """
        if contacts is None:
            return
        contacts = contacts.copy()
        self._update_feature_names(contacts)
        if self._use_thickness:
            contacts['val'] = np.nan
            for k, v in self._stratigraphic_value().items():
                contacts.loc[contacts['name'] == k,'val'] = v

            self._contacts = contacts.loc[~np.isnan(contacts['val']),['X','Y','Z','feature_name','val']]
        if not self._use_thickness:
            contacts['interface'] = np.nan
            for sg in self._stratigraphic_order:
                interface_val = 0
                for g in reversed(sg):
                    contacts.loc[contacts['name'] == g,'interface'] = interface_val
                    interface_val+=1
            self._contacts = contacts.loc[~np.isnan(contacts['interface']),['X','Y','Z','feature_name','interface']]

    
    @property
    def contact_orientations(self):
        contact_orientations = self._contact_orientations.copy()
        #scale
        contact_orientations.loc[~np.isnan(contact_orientations['nz']),['nx', 'ny', 'nz']]*=self.vector_scale*\
                                                                contact_orientations['polarity'].to_numpy()[:,None]
        if self._gradient:
            contact_orientations.rename(columns={'nx':'gx','ny':'gy','nz':'gz'},inplace=True)
        return contact_orientations

    @contact_orientations.setter
    def contact_orientations(self, contact_orientations):
        """process orientation data so strike and dip are converted to normal vectors
        scale vectors as required and convert polarity

        Returns
        -------
        DataFrame
            DataFrame containing orientation data
        """
        if contact_orientations is None:
            return
        if 'polarity' not in contact_orientations.columns:
            contact_orientations['polarity'] = 1.
        contact_orientations = contact_orientations.copy()
        
        contact_orientations['strike'] = contact_orientations['azimuth'] - 90
        contact_orientations['nx'] = np.nan
        contact_orientations['ny'] = np.nan
        contact_orientations['nz'] = np.nan
        contact_orientations[['nx', 'ny', 'nz']] = strike_dip_vector(contact_orientations['strike'],
                                                                    contact_orientations['dip'])
        self._update_feature_names(contact_orientations)
        self._contact_orientations = contact_orientations[['X','Y','Z','nx','ny','nz','feature_name','polarity']]
    
    @property
    def fault_orientations(self):
        return self._fault_orientations

    @fault_orientations.setter
    def fault_orientations(self, fault_orientations):
        """convert fault orientation data to vectors

        Returns
        -------
        DataFrame
            data frame with x,y,z,gx,gy,gz,coord,feature_name
        """
        if fault_orientations is None:
            return
        fault_orientations = fault_orientations.copy()
        fault_orientations['coord'] = 0
        fault_orientations['gx'] = np.nan
        fault_orientations['gy'] = np.nan
        fault_orientations['gz'] = np.nan
        fault_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(fault_orientations['strike'], fault_orientations['dip'])
        fault_orientations['feature_name'] = fault_orientations['fault_name']
        self._fault_orientations = fault_orientations[['X','Y','Z','gx','gy','gz','coord','feature_name']]

    @property
    def fault_locations(self):
        return self._fault_locations

    @fault_locations.setter
    def fault_locations(self, fault_locations):
        """Convert fault traces into input

        Returns
        -------
        Dataframe
            daaframe with x,y,z,val,coord,feature_name
        """
        if fault_locations is None:
            return
        fault_locations = fault_locations.copy()
        fault_locations['coord'] = 0
        fault_locations['val'] = 0
        fault_locations['feature_name'] = fault_locations['fault_name']
        self._fault_locations = fault_locations[['X','Y','Z','val','feature_name','coord']]
    