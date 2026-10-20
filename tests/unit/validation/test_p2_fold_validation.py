"""Check that the validation gate cannot endorse an unfair comparison."""

import pytest

from scripts.validation.summarise_p2_fold import decision


def rows():
    return [
        {
            "target": target,
            "engine": engine,
            "rule": "gauss4",
            "buffer": 0.0,
            "unknowns": unknowns,
            "median_angle_deg": angle,
            "expected_normal_sha256": target,
        }
        for target in ("S2", "S0")
        for engine, unknowns, angle in (("DFI", 12000, 8.0), ("P2Fold", 9000, 5.0))
    ]


def test_gate_requires_improvement_for_both_targets():
    results = rows()
    assert decision(results)["passed"]
    results[-1]["median_angle_deg"] = 9.0
    assert not decision(results)["passed"]


def test_gate_rejects_larger_p2_mesh():
    results = rows()
    results[1]["unknowns"] = 85000
    with pytest.raises(ValueError, match="nodal budget"):
        decision(results)


def test_gate_rejects_different_fold_geometry():
    results = rows()
    results[1]["expected_normal_sha256"] = "different profile"
    with pytest.raises(ValueError, match="different fold events"):
        decision(results)


def test_gate_rejects_different_buffers():
    results = rows()
    results[1]["buffer"] = 0.2
    with pytest.raises(ValueError, match="different buffers"):
        decision(results)
