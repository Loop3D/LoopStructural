import pytest

from LoopStructural.modelling.core.fault_topology import FaultRelationshipType, FaultTopology
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn


@pytest.fixture
def topology():
    sc = StratigraphicColumn()
    topo = FaultTopology(sc)
    topo.add_fault("f1")
    topo.add_fault("f2")
    topo.add_fault("f3")
    return topo


def test_add_and_remove_fault(topology):
    assert topology.get_faults() == ["f1", "f2", "f3"]
    topology.remove_fault("f2")
    assert topology.get_faults() == ["f1", "f3"]


def test_remove_nonexistent_fault_raises(topology):
    with pytest.raises(ValueError):
        topology.remove_fault("does_not_exist")


def test_add_abutting_relationship(topology):
    topology.add_abutting_relationship("f1", "f2")
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.ABUTTING


def test_add_faulted_relationship(topology):
    topology.add_faulted_relationship("f1", "f2")
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.FAULTED


def test_missing_relationship_returns_none_type(topology):
    assert topology.get_fault_relationship("f1", "f3") == FaultRelationshipType.NONE


def test_relationship_requires_both_faults_registered(topology):
    with pytest.raises(ValueError):
        topology.add_abutting_relationship("f1", "does_not_exist")


def test_get_fault_relationships_does_not_crash_on_unpacking(topology):
    """Regression test: add_abutting/faulted_relationship used to also insert a
    dead string-keyed entry (self.adjacency[fault_name] = []) alongside the real
    tuple-keyed one, which broke the `for (f1, f2), relationship_type in
    self.adjacency.items()` unpacking in get_fault_relationships/get_matrix.
    """
    topology.add_abutting_relationship("f1", "f2")
    topology.add_faulted_relationship("f2", "f3")

    rels_f1 = topology.get_fault_relationships("f1")
    rels_f2 = topology.get_fault_relationships("f2")
    assert rels_f1 == [("f1", "f2", FaultRelationshipType.ABUTTING)]
    assert ("f2", "f3", FaultRelationshipType.FAULTED) in rels_f2


def test_get_matrix(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.add_faulted_relationship("f2", "f3")
    matrix = topology.get_matrix()
    assert matrix.shape == (3, 3)
    assert matrix[0, 1] == 1  # abutting
    assert matrix[1, 2] == 2  # faulted
    assert matrix[0, 2] == 0


def test_remove_fault_relationship(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.remove_fault_relationship("f1", "f2")
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.NONE


def test_remove_fault_relationship_reversed_order(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.remove_fault_relationship("f2", "f1")
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.NONE


def test_remove_nonexistent_relationship_raises(topology):
    with pytest.raises(ValueError):
        topology.remove_fault_relationship("f1", "f2")


def test_change_relationship_type(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.change_relationship_type("f1", "f2", FaultRelationshipType.FAULTED)
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.FAULTED


def test_change_relationship_type_requires_existing_relationship(topology):
    with pytest.raises(ValueError):
        topology.change_relationship_type("f1", "f2", FaultRelationshipType.FAULTED)


def test_update_fault_relationship_to_none_removes_it(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.update_fault_relationship("f1", "f2", FaultRelationshipType.NONE)
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.NONE
    assert ("f1", "f2") not in topology.adjacency


def test_remove_fault_clears_its_relationships(topology):
    topology.add_abutting_relationship("f1", "f2")
    topology.add_faulted_relationship("f2", "f3")
    topology.remove_fault("f2")
    assert topology.get_fault_relationship("f1", "f2") == FaultRelationshipType.NONE
    assert topology.get_fault_relationship("f2", "f3") == FaultRelationshipType.NONE


def test_stratigraphy_fault_relationship(topology):
    topology.add_stratigraphy_fault_relationship("unitA", "f1")
    assert topology.get_fault_stratigraphic_relationship("unitA", "f1") is True
    assert topology.get_fault_stratigraphic_relationship("unitB", "f1") is False


def test_update_and_remove_stratigraphy_fault_relationship(topology):
    topology.add_stratigraphy_fault_relationship("unitA", "f1")
    topology.update_fault_stratigraphy_relationship("unitA", "f1", flag=False)
    assert topology.get_fault_stratigraphic_relationship("unitA", "f1") is False

    topology.update_fault_stratigraphy_relationship("unitA", "f1", flag=True)
    assert topology.get_fault_stratigraphic_relationship("unitA", "f1") is True

    topology.remove_fault_stratigraphy_relationship("unitA", "f1")
    assert topology.get_fault_stratigraphic_relationship("unitA", "f1") is False


def test_to_dict_from_dict_round_trip(topology):
    """Regression test: update_from_dict used to iterate adjacency.values()
    instead of .items() and only ever called add_abutting_relationship,
    dropping the FAULTED/ABUTTING distinction and misreading the adjacency
    entries entirely, so to_dict/from_dict was not actually round-trippable.
    """
    topology.add_abutting_relationship("f1", "f2")
    topology.add_faulted_relationship("f2", "f3")
    topology.add_stratigraphy_fault_relationship("unitA", "f1")

    data = topology.to_dict()
    restored = FaultTopology.from_dict(
        {**data, "stratigraphic_column": topology.stratigraphic_column}
    )

    assert restored.get_faults() == topology.get_faults()
    assert restored.adjacency == topology.adjacency
    assert restored.get_fault_relationship("f1", "f2") == FaultRelationshipType.ABUTTING
    assert restored.get_fault_relationship("f2", "f3") == FaultRelationshipType.FAULTED
    assert restored.get_fault_stratigraphic_relationship("unitA", "f1") is True
