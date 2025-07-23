from click import group
from matplotlib.pyplot import cla

from LoopStructural.modelling.core import stratigraphic_column
from ..features.fault import FaultSegment
import enum
import numpy as np
class FaultRelationshipType(enum.Enum):
    ABUTTING = "abutting"
    FAULTED = "faulted"

class FaultTopology:
    """A graph representation of the relationships between faults and the
     relationship with stratigraphic units.
    """
    def __init__(self, stratigraphic_column: 'StratigraphicColumn'):
        self.faults = []
        self.stratigraphic_column = stratigraphic_column
        self.adjacency = {}
        self.stratigraphy_fault_relationships = {}
    def add_fault(self, fault: FaultSegment):
        """
        Adds a fault to the fault topology.
        """
        if not isinstance(fault, FaultSegment):
            raise TypeError("Expected a Fault instance.")
        self.faults.append(fault)

    def add_abutting_relationship(self, fault_name: str, abutting_fault: str):
        """
        Adds an abutting relationship between two faults.
        """
        if fault_name not in self.faults or abutting_fault not in self.faults:
            raise ValueError("Both faults must be part of the fault topology.")

        if fault_name not in self.adjacency:
            self.adjacency[fault_name] = []

        self.adjacency[(fault_name, abutting_fault)] = FaultRelationshipType.ABUTTING
    def add_stratigraphy_fault_relationship(self, unit_name:str, fault_name: str):
        """
        Adds a relationship between a stratigraphic unit and a fault.
        """
        if fault_name not in self.faults:
            raise ValueError("Fault must be part of the fault topology.")

        group = self.stratigraphic_column.get_group_for_unit_name(unit_name)
        if group is None:
            raise ValueError(f"No stratigraphic group found for unit name: {unit_name}")
        if unit_name not in self.stratigraphy_fault_relationships:
            self.stratigraphy_fault_relationships[unit_name] = []
        if fault_name not in self.stratigraphy_fault_relationships[unit_name]:
            self.stratigraphy_fault_relationships[unit_name].append(fault_name)

    def add_faulted_relationship(self, fault_name: str, faulted_fault_name: str):
        """
        Adds a faulted relationship between two faults.
        """
        if fault_name not in self.faults or faulted_fault_name not in self.faults:
            raise ValueError("Both faults must be part of the fault topology.")

        if fault_name not in self.adjacency:
            self.adjacency[fault_name] = []

        self.adjacency[(fault_name, faulted_fault_name)] = FaultRelationshipType.FAULTED
    def remove_fault_relationship(self, fault_name: str, related_fault_name: str):
        """
        Removes a relationship between two faults.
        """
        if (fault_name, related_fault_name) in self.adjacency:
            del self.adjacency[(fault_name, related_fault_name)]
        elif (related_fault_name, fault_name) in self.adjacency:
            del self.adjacency[(related_fault_name, fault_name)]
        else:
            raise ValueError(f"No relationship found between {fault_name} and {related_fault_name}.")
    def change_relationship_type(self, fault_name: str, related_fault_name: str, new_relationship_type: FaultRelationshipType):
        """
        Changes the relationship type between two faults.
        """
        if (fault_name, related_fault_name) in self.adjacency:
            self.adjacency[(fault_name, related_fault_name)] = new_relationship_type
        elif (related_fault_name, fault_name) in self.adjacency:
            self.adjacency[(related_fault_name, fault_name)] = new_relationship_type
        else:
            raise ValueError(f"No relationship found between {fault_name} and {related_fault_name}.")
    
    def get_fault_relationships(self, fault_name: str):
        """
        Returns a list of relationships for a given fault.
        """
        relationships = []
        for (f1, f2), relationship_type in self.adjacency.items():
            if f1 == fault_name or f2 == fault_name:
                relationships.append((f1, f2, relationship_type))
        return relationships
    def get_faults(self):
        """
        Returns a list of all faults in the topology.
        """
        return self.faults
    
    def get_stratigraphy_fault_relationships(self):
        """
        Returns a dictionary of stratigraphic unit to fault relationships.
        """
        return self.stratigraphy_fault_relationships
    def get_fault_stratigraphic_unit_relationships(self):
        units_group_pairs = self.stratigraphic_column.get_units_group_pairs()
        matrix = np.zeros((len(self.faults), len(units_group_pairs)), dtype=int)
        for i, fault in enumerate(self.faults):
            for j, (unit_name, group) in enumerate(units_group_pairs):
                if unit_name in self.stratigraphy_fault_relationships and fault.name in self.stratigraphy_fault_relationships[group]:
                    matrix[i, j] = 1
        return matrix
    def remove_fault_stratigraphy_relationship(self, unit_name: str, fault_name: str):
        """
        Removes a relationship between a stratigraphic unit and a fault.
        """
        group = self.stratigraphic_column.get_group_for_unit_name(unit_name)
        if group is None: 
            raise ValueError(f"No stratigraphic group found for unit name: {unit_name}")
        
        if group in self.stratigraphy_fault_relationships:
            if fault_name in self.stratigraphy_fault_relationships[group]:
                self.stratigraphy_fault_relationships[group].remove(fault_name)
                if not self.stratigraphy_fault_relationships[group]:
                    del self.stratigraphy_fault_relationships[group]
            else:
                raise ValueError(f"Fault {fault_name} not found in relationships for unit {unit_name}.")
        else:
            raise ValueError(f"Unit {unit_name} not found in stratigraphy fault relationships.")
    def get_matrix(self):
        """
        Returns a matrix representation of the fault relationships.
        """
        matrix = np.zeros((len(self.faults), len(self.faults)), dtype=int)
        for (fault_name, related_fault_name), relationship_type in self.adjacency.items():
            fault_index = self.faults.index(next(f for f in self.faults if f == fault_name))
            related_fault_index = self.faults.index(next(f for f in self.faults if f == related_fault_name))
            if relationship_type == FaultRelationshipType.ABUTTING:
                matrix[fault_index, related_fault_index] = 1
            elif relationship_type == FaultRelationshipType.FAULTED:
                matrix[fault_index, related_fault_index] = 2
        return matrix

    def to_dict(self):
        """
        Returns a dictionary representation of the fault topology.
        """
        return {
            "faults": self.faults,
            "adjacency": self.adjacency,
            "stratigraphy_fault_relationships": self.stratigraphy_fault_relationships,
        }
    
    def update_from_dict(self, data):
        """
        Updates the fault topology from a dictionary representation.
        """
        self.faults.extend(data.get("faults", []))
        adjacency = data.get("adjacency", {})
        stratigraphy_fault_relationships = data.get("stratigraphy_fault_relationships", {})
        for (fault,abutting_fault) in adjacency.values():
            if fault not in self.faults:
                self.add_fault(fault)
            if abutting_fault not in self.faults:
                self.add_fault(abutting_fault)
            self.add_abutting_relationship(fault, abutting_fault)
        for unit_name, fault_names in stratigraphy_fault_relationships.items():
            for fault_name in fault_names:
                if fault_name not in self.faults:
                    self.add_fault(fault_name)
                self.add_stratigraphy_fault_relationship(unit_name, fault_name)
    
    @classmethod
    def from_dict(cls, data):
        """
        Creates a FaultTopology instance from a dictionary representation.
        """
        from .stratigraphic_column import StratigraphicColumn
        stratigraphic_column = data.get("stratigraphic_column",None)
        if not isinstance(stratigraphic_column, StratigraphicColumn):
            if isinstance(stratigraphic_column, dict):
                stratigraphic_column = StratigraphicColumn.from_dict(stratigraphic_column)
        elif not isinstance(stratigraphic_column, StratigraphicColumn):
            raise TypeError("Expected 'stratigraphic_column' to be a StratigraphicColumn instance or dict.")
        
        topology = cls(stratigraphic_column)
        topology.update_from_dict(data)
        return topology