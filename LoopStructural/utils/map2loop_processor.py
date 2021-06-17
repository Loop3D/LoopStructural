import pandas as pd
import numpy as np
import os
import logging
import networkx
from scipy.stats import truncnorm

from LoopStructural.utils import strike_dip_vector
from LoopStructural.utils import getLogger
logger = getLogger(__name__)
class ProcessInputData:
    def __init__(   self,  
                    contacts, 
                    orientations, 
                    fault_displacements, 
                    fault_orientations, 
                    fault_locations, 
                    fault_dimensions, 
                    fault_graph,
                    stratigraphic_order,
                    intrusions,
                    fault_stratigraphy
                ):
 
        self._contacts = contacts
        self._orientations = orientations
        self._fault_displacements = fault_displacements
        self._fault_orientations = fault_orientations
        self._fault_locations = fault_locations
        self._fault_dimensions = fault_dimensions
        self._fault_graph = fault_graph
        self._fault_stratigraphy = fault_stratigraphy
        self._stratigraphic_order = stratigraphic_order
        self._intrusions = intrusions
        self._thicknesses = {}
        self._data = None
    
    @property
    def data(self):
        self._update()
        return self._data
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
  

    
    def _update(self):
        dataframes = []
        dataframes.append(self._process_contacts())
        dataframes.append(self._process_orientations())
        dataframes.append(self._process_fault_orientations())
        dataframes.append(self._process_fault_locations())
        self._data = pd.concat(dataframes)
        self._data.reset_index(inplace=False)
    
    def _process_contacts(self):
        value = 0
        stratigraphic_values = {}
        for sg in self._stratigraphic_order:
            for g in reversed(sg):
                stratigraphic_values[g] = value
                value+=self._thicknesses[g]
        use_thickness = flags.get('use_thickness',True)
        if use_thickness:
            contacts['val'] = np.nan
            for o in strat_val:
                contacts.loc[contacts['formation'] == o, 'val'] = strat_val[o]
    if use_thickness == False:
        contacts['interface'] = np.nan
        interface_val = 0
        for u in contacts['formation'].unique():
            contacts.loc[contacts['formation'] == u,'interface'] = interface_val
            interface_val+=1
    tangents['feature_name'] = tangents['group']
    contact_orientations['feature_name'] = None
    contacts['feature_name'] = None
    for g in groups['group'].unique():
        val = 0
        for c in groups.loc[groups['group'] == g, 'code']:
            contact_orientations.loc[contact_orientations['formation'] == c, 'feature_name'] = supergroups[g]
            contacts.loc[contacts['formation'] == c, 'feature_name'] = supergroups[g]
    def _process_orientations(self):
        if self._gradient:
            contact_orientations['strike'] = contact_orientations['azimuth'] - 90
            contact_orientations['gx'] = np.nan
            contact_orientations['gy'] = np.nan
            contact_orientations['gz'] = np.nan
            contact_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(contact_orientations['strike'],
                                                                        contact_orientations['dip'])
                                                                        *self._vector_scale   
        if np.sum(contact_orientations['polarity']==0) >0 and np.sum(contact_orientations['polarity']==-1)==0:
            # contact_orientations['polarity']+=1
            contact_orientations.loc[contact_orientations['polarity']==0]=-1
        if not gradient:
            from LoopStructural.utils.helper import strike_dip_vector
            contact_orientations['strike'] = contact_orientations['azimuth'] - 90
            contact_orientations['nx'] = np.nan
            contact_orientations['ny'] = np.nan
            contact_orientations['nz'] = np.nan
            contact_orientations[['nx', 'ny', 'nz']] = strike_dip_vector(contact_orientations['strike'],
                                                                        contact_orientations['dip'])*vector_scale *contact_orientations['polarity'].to_numpy()[:,None]
        contact_orientations.drop(['strike', 'dip', 'azimuth'], inplace=True, axis=1)
    def _process_fault_orientations(self):

    def _process_fault_locations(self):

    def load_unit_unit_thicknesses(self,thickness_file, thickness_probabilities = False):
        i = 0
        thickness = {}
        max_thickness = 0
    
        with open(thickness_file) as file:
            for l in file:
                if i>=1:
                    linesplit = l.split(',')
                    if thickness_probabilities:
                        std = float(linesplit[2])
                        mean = float(linesplit[1])
                        if  np.isnan(std):
                            std = float(linesplit[1])
                        a = (0-mean) / std
                        b = 100.
                        thickness[linesplit[0]] = (truncnorm.rvs(a,b,size=1)*std)[0]+mean
                    else:
                        thickness[linesplit[0]] = float(linesplit[1])
                    
                    # normalise the thicknesses
                    if float(linesplit[1]) > max_thickness:
                        max_thickness=float(linesplit[1])
        #             print(l.split(',')[1])
                i+=1
        self.thicknesses = thicknesses


    groups = pd.read_csv(m2l_directory + '/tmp/all_sorts_clean.csv', index_col=0)
    contact_orientations = pd.read_csv(m2l_directory + '/output/orientations_clean.csv')
    # formation_thickness = pd.read_csv)
    contacts = pd.read_csv(m2l_directory + '/output/contacts_clean.csv')
    displacements = pd.read_csv(m2l_directory + '/output/fault_displacements3.csv')
    fault_orientations = pd.read_csv(m2l_directory + '/output/fault_orientations.csv')
    fault_locations = pd.read_csv(m2l_directory + '/output/faults.csv')
    fault_fault_relations = pd.read_csv(m2l_directory + '/output/fault-fault-relationships.csv')
    fault_strat_relations = pd.read_csv(m2l_directory + '/output/group-fault-relationships.csv')
    fault_dimensions = pd.read_csv(m2l_directory + '/output/fault_dimensions.csv')
    fault_graph = networkx.read_gml(m2l_directory + '/tmp/fault_network.gml')