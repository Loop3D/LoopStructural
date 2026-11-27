"""The 94-fault benchmark model of spec 007: its fault layout is fixed (no model build)."""

import numpy as np

from benchmarks.fault_evaluation import many_faults


def test_fault_layout_is_deterministic():
    specs = many_faults.fault_specs()
    assert len(specs) == 94
    assert sum(s["kind"] == "cutting" for s in specs) == 40
    # A fault is cut only by younger (lower index) cutting faults, at most 3 of them.
    for k, s in enumerate(specs):
        assert all(m < k and specs[m]["kind"] == "cutting" for m in s["cut_by"])
        assert len(s["cut_by"]) <= many_faults.MAX_CUTTERS
    assert [s["cut_by"] for s in many_faults.fault_specs()] == [s["cut_by"] for s in specs]


def test_compare_reports_changed_and_missing_arrays(tmp_path):
    path = tmp_path / "ref.npz"
    np.savez(path, a=np.array([1.0, np.nan]), b=np.array([1]))
    assert many_faults.compare({"a": np.array([1.0, np.nan]), "b": np.array([1])}, path) == []
    assert many_faults.compare({"a": np.array([1.0, 2.0])}, path) == ["a", "b"]


def test_fault_summary_is_none_when_no_fault_call_is_timed():
    summary = many_faults._fault_summary({"isolated": [], "cutting": []})
    assert summary == {"isolated_mean_ms": None, "cutting_mean_ms": None, "cutting_max_ms": None}


def test_fault_summary_of_timed_calls():
    summary = many_faults._fault_summary({"isolated": [1.0, 3.0], "cutting": [2.0, 6.0]})
    assert summary == {"isolated_mean_ms": 2.0, "cutting_mean_ms": 4.0, "cutting_max_ms": 6.0}
