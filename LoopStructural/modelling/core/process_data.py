import pandas as pd
import numpy as np
import os
import logging
import networkx
from scipy.stats import truncnorm
import copy
from LoopStructural.utils import strike_dip_vector
from LoopStructural.utils import getLogger
logger = getLogger(__name__)
class ProcessInputData:
    def __init__(   self, 
                    contacts, 
                    contact_orientations, 
                    stratigraphic_order,
                    fault_orientations = None, 
                    fault_locations = None, 
                    fault_dimensions = None, 
                    fault_graph = None,
                    intrusions = None,
                    fault_stratigraphy = None,
                    thicknesses = None,
                    colours = None
                ):
        """

        Parameters
        ----------
        contacts : DataFrame
            x,y,z,formation
        contact_orientations : DataFrame
            x,y,z,strike,dip, formation
        stratigraphic_order : nested list
            a nested list e.g. [['a','b','c'],['d','e']]
            a->b->c are the youngest supergroup, d->e older.
        fault_orientations : DataFrame, optional
            data frame with x,y,z,strike,fault_name, by default None
        fault_locations : DataFrame, optional
            data frame with x,y,z,fault_name, by default None
        fault_dimensions : DataFrame, optional
            dataframe with fault_name, InfluenceDistance,SlipDistance,ExtentDistance,
             Displacement, colour, X,Y,Z by default None
        fault_graph : graph, optional
            graph describing fault fault relationship, by default None
        intrusions : list, optional
            list containing any intrusions, assumes data in contacts/orientations, by default None
        fault_stratigraphy : graph, optional
            fault stratigraphy intersection, by default None
        thicknesses : dict, optional
            dictionary with stratigraphic thicknesses, by default None
        colours : dict, optional
            dictionary with stratigraphic unit colours, by default None

        Notes
        ------
        
        """
        
        self._stratigraphic_order = stratigraphic_order
        self._thicknesses = thicknesses
        if self.thicknesses is None:
            self._use_thickness = False
        else:
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
        self._fault_dimensions = None
        self.fault_dimensions = fault_dimensions
        self._fault_graph = None
        self.fault_graph = fault_graph
        self._fault_stratigraphy = fault_stratigraphy
        self._intrusions = intrusions
        self._thicknesses = thicknesses
        self._data = None
        self._colours = {}
        self.colours = colours
        # flags
 

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
        for i, sg in enumerate(self._stratigraphic_order):
            stratigraphic_column['supergroup_{}'.format(i)] = {}
            for g in reversed(sg):
                stratigraphic_column['supergroup_{}'.format(i)][g] = {'max': val[g]+self.thicknesses[g], 'min': val[g] , 'id': unit_id, 'colour':self.colours[g]}
        # add faults into the column 
        stratigraphic_column['faults'] = self.fault_dimensions
        return stratigraphic_column

    @property
    def fault_dimensions(self):
        return self._fault_dimensions

    @fault_dimensions.setter
    def fault_dimensions(self, fault_dimensions):
        if fault_dimensions is None:
            return
        fault_dimensions = fault_dimensions.copy()
        if 'X' not in fault_dimensions.columns or \
           'Y' not in fault_dimensions.columns or \
           'Z' not in fault_dimensions:
            fault_dimensions['X'] = np.nan
            fault_dimensions['Y'] = np.nan
            fault_dimensions['Z'] = np.nan
            for fname in fault_dimensions.index:
                pts = self.fault_locations.loc[self.fault_locations['feature_name'] == fname,
                                                        ['X','Y','Z']]
                fault_dimensions.loc[fname,['X','Y','Z']] = np.mean(pts,axis=0)
        if 'nx' not in fault_dimensions.columns or \
           'ny' not in fault_dimensions.columns or \
           'nz' not in fault_dimensions:
            fault_dimensions['nx'] = np.nan
            fault_dimensions['ny'] = np.nan
            fault_dimensions['nz'] = np.nan

            for fname in fault_dimensions.index:
                pts = self.fault_orientations.loc[self.fault_orientations['feature_name'] == fname,
                                                        ['gx','gy','gz']]
                fault_dimensions.loc[fname,['nx','ny','nz']] = np.mean(pts,axis=0)
        if 'sx' not in fault_dimensions.columns or \
           'sy' not in fault_dimensions.columns or \
           'sz' not in fault_dimensions:
            fault_dimensions['sx'] = np.nan
            fault_dimensions['sy'] = np.nan
            fault_dimensions['sz'] = np.nan
            for fname in fault_dimensions.index:
                if  'sx' not in self.fault_orientations.columns or \
                    'sy' not in self.fault_orientations.columns or \
                    'sz' not in self.fault_orientations:
                    # if we don't know slip assum its down
                    fault_dimensions.loc[fname,['sx','sy','sz']] = np.array([0.,0.,-1.])
                else:
                    pts = self.fault_orientations.loc[self.fault_orientations['feature_name'] == fname,
                                                            ['sx','sy','sz']]
                    fault_dimensions.loc[fname,['sx','sy','sz']] = np.mean(pts,axis=0)

        self._fault_dimensions = fault_dimensions

    @property
    def fault_graph(self):
        return self._fault_graph
    def draw_fault_graph(self):
        try:
            import networkx
        
            pos=networkx.nx_agraph.graphviz_layout(self.fault_graph,prog="dot")
            networkx.draw(self.fault_graph,pos,with_labels=True)
        except ImportError:
            logger.error("dependencies missing can't plot graph")
    
    @fault_graph.setter
    def fault_graph(self, fault_graph):
        fault_graph = copy.copy(fault_graph)
        # add attributes from the fault_dimensions table
        for n in fault_graph.nodes:
            fault_graph.nodes[n]['fault_influence'] = self.fault_dimensions.loc[n,'InfluenceDistance']
            fault_graph.nodes[n]['fault_vectical_radius'] = self.fault_dimensions.loc[n,'SlipDistance']
            fault_graph.nodes[n]['fault_extent'] = self.fault_dimensions.loc[n,'ExtentDistance']
            fault_graph.nodes[n]['colour'] = self.fault_dimensions.loc[n,'colour']
            fault_graph.nodes[n]['fault_normal'] = self.fault_dimensions.loc[n,['nx','ny','nz']].to_numpy().astype(float)
            fault_graph.nodes[n]['fault_slip_vector'] = self.fault_dimensions.loc[n,['sx','sy','sz']].to_numpy().astype(float)
            fault_graph.nodes[n]['displacement'] = self.fault_dimensions.loc[n,'displacement']
            fault_graph.nodes[n]['fault_center'] = self.fault_dimensions.loc[n,['X','Y','Z']].to_numpy().astype(float)
            fault_graph.nodes[n]['faultfunction'] = 'BaseFault'
        self._fault_graph = fault_graph

    @property
    def fault_names(self):
        return list(self.fault_dimensions['fault_name'])
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
            dataframes.append(self.fault_orientations)
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
        for sg in self._stratigraphic_order:
            value = 0 #reset for each supergroup
            for g in reversed(sg):
                stratigraphic_value[g] = value
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
            for i, sg in enumerate(self._stratigraphic_order):
                for g in sg:
                    dataframe.loc[dataframe['formation']==g,'feature_name'] = 'supergoup_{}'.format(i)
    
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
                contacts.loc[contacts['formation'] == k,'val'] = v

            self._contacts = contacts.loc[~np.isnan(contacts['val']),['X','Y','Z','feature_name','val']]
        if not self._use_thickness:
            contacts['interface'] = np.nan
            for sg in self._stratigraphic_order:
                interface_val = 0
                for g in reversed(sg):
                    contacts.loc[contacts['formation'] == g,'interface'] = interface_val
                    interface_val+=1
            self._contacts = contacts.loc[~np.isnan(contacts['interface']),['X','Y','Z','feature_name','interface']]

    
    @property
    def contact_orientations(self):
        return self._contact_orientations

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
        contact_orientations = contact_orientations.copy()
        contact_orientations['strike'] = contact_orientations['azimuth'] - 90
        contact_orientations['nx'] = np.nan
        contact_orientations['ny'] = np.nan
        contact_orientations['nz'] = np.nan
        contact_orientations[['nx', 'ny', 'nz']] = strike_dip_vector(contact_orientations['strike'],
                                                                    contact_orientations['dip'])
        self._update_feature_names(contact_orientations)
        contact_orientations.loc[~np.isnan(contact_orientations['nz']),['nx', 'ny', 'nz']]*=self._vector_scale*\
                                                                contact_orientations['polarity'].to_numpy()[:,None]
        contact_orientations.drop(['strike', 'dip', 'azimuth'], inplace=True, axis=1)
        if self._gradient:
            contact_orientations.rename(columns={'nx':'gx','ny':'gy','nz':'gz'},inplace=True)
        self._contact_orientations = contact_orientations[['X','Y','Z','nx','ny','nz','feature_name']]
    
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
    