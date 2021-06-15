import pandas as pd
import numpy as np
import os
import logging
import networkx
from scipy.stats import truncnorm

from LoopStructural.utils import getLogger
logger = getLogger(__name__)
class FromMap2Loop:
    def __init__(self, contacts, displacements, fault_orientations, fault_locations, fault_dimensions, fault_graph ):
        self._contacts = contacts
        self._displacements = displacements
        self._fault_orientations = fault_orientations
        self._fault_locations = fault_locations
        self._fault_dimensions = fault_dimensions
        self._fault_graph = fault_graph
        self._thicknesses = {}
        self._data = pd.DataFrame()
    @property
    def thicknesses(self):
        return self._thicknesses
    
    @thicknesses.setter
    def thicknesses(self,thicknesses):
        self._thicknesses = thicknesses

    def calculate_unit_thicknesses(self,thickness_file, thickness_probabilities = False):
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