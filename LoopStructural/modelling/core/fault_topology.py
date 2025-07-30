from ..features.fault import FaultSegment
from ...utils import Observable
from .stratigraphic_column import StratigraphicColumn
import enum
import numpy as np
class FaultRelationshipType(enum.Enum):
    ABUTTING = "abutting"
    FAULTED = "faulted"
    NONE = "none"

class FaultTopology(Observable['FaultTopology']):
    """A graph representation of the relationships between faults and the
     relationship with stratigraphic units.
    """
    def __init__(self, stratigraphic_column: 'StratigraphicColumn'):
        super().__init__()
        self.faults = []
        self.stratigraphic_column = stratigraphic_column
        self.adjacency = {}
        self.stratigraphy_fault_relationships = {}
    def add_fault(self, fault: FaultSegment):
        """
        Adds a fault to the fault topology.
        """
        if not isinstance(fault, str):
            raise TypeError("Expected a fault name.")

        self.faults.append(fault)
        self.notify('fault_added', fault=fault)

    def remove_fault(self, fault: str):
        """
        Removes a fault from the fault topology.
        """
        if fault not in self.faults:
            raise ValueError(f"Fault {fault} not found in the topology.")
        
        self.faults.remove(fault)
        # Remove any relationships involving this fault
        self.adjacency = {k: v for k, v in self.adjacency.items() if fault not in k}
        self.stratigraphy_fault_relationships = {
            k: v for k, v in self.stratigraphy_fault_relationships.items() if k[1] != fault
        }
        self.notify('fault_removed', fault=fault)
        
    def add_abutting_relationship(self, fault_name: str, abutting_fault: str):
        """
        Adds an abutting relationship between two faults.
        """
        if fault_name not in self.faults or abutting_fault not in self.faults:
            raise ValueError("Both faults must be part of the fault topology.")

        if fault_name not in self.adjacency:
            self.adjacency[fault_name] = []

        self.adjacency[(fault_name, abutting_fault)] = FaultRelationshipType.ABUTTING
        self.notify('abutting_relationship_added', {'fault': fault_name, 'abutting_fault': abutting_fault})
    def add_stratigraphy_fault_relationship(self, unit_name:str, fault_name: str):
        """
        Adds a relationship between a stratigraphic unit and a fault.
        """
        if fault_name not in self.faults:
            raise ValueError("Fault must be part of the fault topology.")

        if unit_name is None:
            raise ValueError(f"No stratigraphic group found for unit name: {unit_name}")
        self.stratigraphy_fault_relationships[(unit_name,fault_name)] = True

        self.notify('stratigraphy_fault_relationship_added', {'unit': unit_name, 'fault': fault_name})
    def add_faulted_relationship(self, fault_name: str, faulted_fault_name: str):
        """
        Adds a faulted relationship between two faults.
        """
        if fault_name not in self.faults or faulted_fault_name not in self.faults:
            raise ValueError("Both faults must be part of the fault topology.")

        if fault_name not in self.adjacency:
            self.adjacency[fault_name] = []

        self.adjacency[(fault_name, faulted_fault_name)] = FaultRelationshipType.FAULTED
        self.notify('faulted_relationship_added', {'fault': fault_name, 'faulted_fault': faulted_fault_name})
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
        self.notify('fault_relationship_removed', {'fault': fault_name, 'related_fault': related_fault_name})
    def update_fault_relationship(self, fault_name: str, related_fault_name: str, new_relationship_type: FaultRelationshipType):
        if new_relationship_type == FaultRelationshipType.NONE:
            self.adjacency.pop((fault_name, related_fault_name), None)
        else:
            self.adjacency[(fault_name, related_fault_name)] = new_relationship_type
        self.notify('fault_relationship_updated', {'fault': fault_name, 'related_fault': related_fault_name, 'new_relationship_type': new_relationship_type})
    def change_relationship_type(self, fault_name: str, related_fault_name: str, new_relationship_type: FaultRelationshipType):
        """
        Changes the relationship type between two faults.
        """
        if (fault_name, related_fault_name) in self.adjacency:
            self.adjacency[(fault_name, related_fault_name)] = new_relationship_type

        else:
            raise ValueError(f"No relationship found between {fault_name} and {related_fault_name}.")
        self.notify('relationship_type_changed', {'fault': fault_name, 'related_fault': related_fault_name, 'new_relationship_type': new_relationship_type})
    def get_fault_relationships(self, fault_name: str):
        """
        Returns a list of relationships for a given fault.
        """
        relationships = []
        for (f1, f2), relationship_type in self.adjacency.items():
            if f1 == fault_name or f2 == fault_name:
                relationships.append((f1, f2, relationship_type))
        return relationships
    def get_fault_relationship(self, fault_name: str, related_fault_name: str):
        """
        Returns the relationship type between two faults.
        """
        return self.adjacency.get((fault_name, related_fault_name), FaultRelationshipType.NONE)
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
        units_group_pairs = self.stratigraphic_column.get_group_unit_pairs() 
        matrix = np.zeros((len(self.faults), len(units_group_pairs)), dtype=int)
        for i, fault in enumerate(self.faults):
            for j, (unit_name, _group) in enumerate(units_group_pairs):
                if (unit_name, fault) in self.stratigraphy_fault_relationships:
                    matrix[i, j] = 1

        return matrix
    def get_fault_stratigraphic_relationship(self, unit_name: str, fault:str) -> bool:
        """
        Returns a dictionary of fault to stratigraphic unit relationships.
        """
        if unit_name is None:
            raise ValueError(f"No stratigraphic group found for unit name: {unit_name}")
        if (unit_name, fault) not in self.stratigraphy_fault_relationships:
            return False
        return self.stratigraphy_fault_relationships[(unit_name, fault)]

    def update_fault_stratigraphy_relationship(self, unit_name: str, fault_name: str, flag: bool = True):
        """
        Updates the relationship between a stratigraphic unit and a fault.
        """
        if not flag:
            if (unit_name, fault_name) in self.stratigraphy_fault_relationships:
                del self.stratigraphy_fault_relationships[(unit_name, fault_name)]
        else:
            self.stratigraphy_fault_relationships[(unit_name, fault_name)] = flag

        self.notify('stratigraphy_fault_relationship_updated', {'unit': unit_name, 'fault': fault_name})

    def remove_fault_stratigraphy_relationship(self, unit_name: str, fault_name: str):
        """
        Removes a relationship between a stratigraphic unit and a fault.
        """
        if (unit_name, fault_name) not in self.stratigraphy_fault_relationships:
            raise ValueError(f"No relationship found between unit {unit_name} and fault {fault_name}.")
        else:
            self.stratigraphy_fault_relationships.pop((unit_name, fault_name), None)

        self.notify('stratigraphy_fault_relationship_removed', {'unit': unit_name, 'fault': fault_name})
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
        with self.freeze_notifications():
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
