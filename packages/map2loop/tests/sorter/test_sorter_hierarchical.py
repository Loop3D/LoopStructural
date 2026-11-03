import pandas
import pytest

from map2loop.sorter import (
    SorterAgeBased,
    SorterAlpha,
    SorterHierarchical,
    SorterUseNetworkX,
    relabel_contacts,
    relabel_unit_relationships,
)


def make_units(rows):
    units = pandas.DataFrame(rows, columns=["name", "minAge", "maxAge", "group", "supergroup"])
    units["layerId"] = units.index
    return units


def assert_contiguous(order, units, column):
    """Check that the units of each value of column are next to each other in order"""
    values = units.set_index("name")[column]
    sequence = [values[name] for name in order]
    seen = []
    for value in sequence:
        if seen and seen[-1] == value:
            continue
        assert value not in seen, f"{column} {value} is split in {order}"
        seen.append(value)


def test_relabel_contacts_adds_lengths_and_removes_internal_contacts():
    contacts = pandas.DataFrame(
        {
            "UNITNAME_1": ["A", "B", "C", "A", "X"],
            "UNITNAME_2": ["B", "C", "A", "C", "A"],
            "length": [10.0, 5.0, 2.0, 3.0, 100.0],
        }
    )
    labels = {"A": "G1", "B": "G1", "C": "G2"}
    result = relabel_contacts(contacts, labels)
    assert len(result) == 1
    row = result.iloc[0]
    assert {row["UNITNAME_1"], row["UNITNAME_2"]} == {"G1", "G2"}
    # B-C, C-A and A-C; A-B is in G1 and X is not in labels
    assert row["length"] == pytest.approx(10.0)


def test_relabel_unit_relationships_keeps_direction():
    relationships = pandas.DataFrame(
        {"UNITNAME_1": ["A", "B", "A"], "UNITNAME_2": ["B", "C", "C"]}
    )
    labels = {"A": "G1", "B": "G1", "C": "G2"}
    result = relabel_unit_relationships(relationships, labels)
    assert result.values.tolist() == [["G1", "G2"]]


def test_age_based_keeps_groups_together():
    units = make_units(
        [
            ["A", 1.0, 2.0, "G1", "S1"],
            ["B", 3.0, 4.0, "G1", "S1"],
            ["C", 2.5, 2.6, "G2", "S1"],
        ]
    )
    # With no hierarchy, C is between A and B
    assert SorterAgeBased().sort(units.drop(columns=["group"])) == ["A", "C", "B"]
    order = SorterHierarchical(sorter=SorterAgeBased()).sort(units)
    assert order == ["A", "B", "C"]


def test_two_supergroups_have_separate_orders():
    units = make_units(
        [
            ["A", 0, 0, "", "S1"],
            ["B", 0, 0, "", "S1"],
            ["C", 0, 0, "", "S1"],
            ["D", 0, 0, "", "S2"],
            ["E", 0, 0, "", "S2"],
        ]
    )
    relationships = pandas.DataFrame(
        {
            "UNITNAME_1": ["A", "B", "D", "B"],
            "UNITNAME_2": ["B", "C", "E", "D"],
        }
    )
    sorter = SorterHierarchical(sorter=SorterUseNetworkX(unit_relationships=relationships))
    order = sorter.sort(units)
    assert order == ["A", "B", "C", "D", "E"]
    # the original sorter data is not changed
    assert len(sorter.sorter.unit_relationships) == 4


def test_groups_in_supergroups_with_contacts():
    units = make_units(
        [
            ["A", 0, 0, "G1", "S1"],
            ["B", 0, 0, "G2", "S1"],
            ["C", 0, 0, "G1", "S1"],
            ["D", 0, 0, "G2", "S1"],
            ["E", 0, 0, "G3", "S2"],
            ["F", 0, 0, "G3", "S2"],
        ]
    )
    contacts = pandas.DataFrame(
        {
            "UNITNAME_1": ["A", "C", "B", "B", "D", "E", "A"],
            "UNITNAME_2": ["C", "B", "D", "E", "F", "F", "B"],
            "length": [100.0, 50.0, 80.0, 500.0, 30.0, 90.0, 400.0],
        }
    )
    order = SorterHierarchical(sorter=SorterAlpha(contacts=contacts)).sort(units)
    assert sorted(order) == sorted(units["name"])
    assert_contiguous(order, units, "supergroup")
    assert_contiguous(order, units, "group")


def test_unit_with_no_group_is_one_item():
    units = make_units(
        [
            ["A", 1.0, 2.0, "G1", ""],
            ["B", 3.0, 4.0, "G1", ""],
            ["C", 2.5, 2.6, None, ""],
            ["D", 5.0, 6.0, "G2", ""],
        ]
    )
    order = SorterHierarchical(sorter=SorterAgeBased()).sort(units)
    assert order == ["A", "B", "C", "D"]


def test_group_with_no_supergroup_is_one_item():
    units = make_units(
        [
            ["A", 1.0, 2.0, "G1", "S1"],
            ["B", 9.0, 10.0, "G2", "S1"],
            ["C", 4.0, 5.0, "G3", ""],
            ["D", 5.0, 6.0, "G3", ""],
            ["E", 3.0, 4.0, "G4", "S2"],
        ]
    )
    order = SorterHierarchical(sorter=SorterAgeBased()).sort(units)
    # S1 spans 1 to 10 (mean 5.5), G3 spans 4 to 6 (mean 5), S2 spans 3 to 4
    assert order == ["E", "C", "D", "A", "B"]


def test_sorter_failure_keeps_all_units():
    units = make_units(
        [
            ["A", 0, 0, "G1", ""],
            ["B", 0, 0, "G1", ""],
            ["C", 0, 0, "G2", ""],
        ]
    )
    # no contacts in G1, so SorterAlpha can not sort A and B
    contacts = pandas.DataFrame(
        {"UNITNAME_1": ["B"], "UNITNAME_2": ["C"], "length": [10.0]}
    )
    order = SorterHierarchical(sorter=SorterAlpha(contacts=contacts)).sort(units)
    assert sorted(order) == ["A", "B", "C"]
    assert_contiguous(order, units, "group")


def test_postprocess_is_applied_to_each_step():
    units = make_units(
        [
            ["A", 1.0, 2.0, "G1", ""],
            ["B", 3.0, 4.0, "G1", ""],
            ["C", 5.0, 6.0, "G2", ""],
            ["D", 7.0, 8.0, "G2", ""],
        ]
    )
    calls = []

    def reverse(order, labels):
        calls.append(sorted(set(labels.values())))
        return list(reversed(order))

    order = SorterHierarchical(sorter=SorterAgeBased(), postprocess=reverse).sort(units)
    assert order == ["D", "C", "B", "A"]
    assert ["G1", "G2"] in calls


def failing_step_units_and_contacts():
    units = make_units(
        [
            ["A", 0, 0, "G1", ""],
            ["B", 0, 0, "G1", ""],
            ["C", 0, 0, "G2", ""],
        ]
    )
    # no contacts in G1, so SorterAlpha can not sort A and B
    contacts = pandas.DataFrame({"UNITNAME_1": ["B"], "UNITNAME_2": ["C"], "length": [10.0]})
    return units, contacts


def test_strict_raises_the_error_of_the_inner_sorter():
    units, contacts = failing_step_units_and_contacts()
    sorter = SorterHierarchical(sorter=SorterAlpha(contacts=contacts), strict=True)
    with pytest.raises(ValueError, match="contacts"):
        sorter.sort(units)


def test_not_strict_is_the_default_and_keeps_all_units():
    units, contacts = failing_step_units_and_contacts()
    sorter = SorterHierarchical(sorter=SorterAlpha(contacts=contacts))
    assert sorter.strict is False
    assert sorted(sorter.sort(units)) == ["A", "B", "C"]


def test_errors_that_are_not_sorting_errors_are_raised_even_when_not_strict():
    class Broken(SorterAgeBased):
        def sort(self, units):
            raise AttributeError("a bug in the sorter")

    units = make_units([["A", 1.0, 2.0, "G1", ""], ["B", 3.0, 4.0, "G2", ""]])
    with pytest.raises(AttributeError):
        SorterHierarchical(sorter=Broken()).sort(units)


def test_nested_hierarchical_sorter_is_a_type_error():
    inner = SorterHierarchical(sorter=SorterAgeBased())
    with pytest.raises(TypeError, match="SorterHierarchical"):
        SorterHierarchical(sorter=inner)


def test_sort_checks_the_type_of_units():
    with pytest.raises(Exception, match="(?i)units|DataFrame"):
        SorterHierarchical(sorter=SorterAgeBased()).sort(["A", "B"])


def test_single_unit():
    units = make_units([["A", 1.0, 2.0, "G1", "S1"]])
    assert SorterHierarchical(sorter=SorterAgeBased()).sort(units) == ["A"]
    assert SorterHierarchical(sorter=SorterAgeBased(), strict=True).sort(units) == ["A"]


def test_empty_units():
    units = make_units([])
    assert SorterHierarchical(sorter=SorterAgeBased()).sort(units) == []


def test_units_without_the_name_column_is_a_value_error():
    with pytest.raises(ValueError, match="name"):
        SorterHierarchical(sorter=SorterAgeBased()).sort(pandas.DataFrame())
