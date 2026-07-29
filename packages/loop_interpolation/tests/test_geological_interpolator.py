import numpy as np
import pytest

from loop_common.interfaces.representation import BaseRepresentation
from loop_interpolation import GeologicalInterpolator
from loop_interpolation.constraints import ValueConstraint, GradientConstraint


def test_get_data_locations(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    locations = interpolator.get_data_locations()
    assert np.sum(locations - data[["X", "Y", "Z"]].to_numpy()) == 0


def test_get_value_constraints(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    val = interpolator.get_value_constraints()
    assert np.sum(val - data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()) == 0


def test_get_norm_constraints(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    val = interpolator.get_norm_constraints()
    assert (
        np.sum(
            val - data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
        )
        == 0
    )


def test_reset(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    interpolator.clean()
    assert interpolator.get_data_locations().shape[0] == 0
    assert not interpolator.up_to_date


def test_interpolator_is_base_representation(interpolator):
    assert isinstance(interpolator, BaseRepresentation)


def test_geological_interpolator_from_dict_delegates_to_factory(monkeypatch):
    from loop_interpolation._interpolator_factory import InterpolatorFactory

    payload = {"type": "FDI", "custom": "value"}
    sentinel = object()
    calls = []

    def _fake_from_dict(data):
        calls.append(data)
        return sentinel

    monkeypatch.setattr(InterpolatorFactory, "from_dict", _fake_from_dict)

    result = GeologicalInterpolator.from_dict(payload)

    assert result is sentinel
    assert calls == [payload]


class _MinimalGeologicalInterpolator(GeologicalInterpolator):
    def __init__(self):
        super().__init__()

    def set_nelements(self, nelements: int) -> int:
        return nelements

    @property
    def n_elements(self) -> int:
        return 0

    def set_region(self, **kwargs):
        return None

    def setup_interpolator(self, **kwargs):
        return None

    def solve_system(self, solver, solver_kwargs: dict = {}) -> bool:
        return True

    def update(self) -> bool:
        return True

    def evaluate_value(self, locations: np.ndarray):
        return np.zeros(np.asarray(locations).shape[0])

    def evaluate_gradient(self, locations: np.ndarray):
        locations = np.asarray(locations)
        return np.zeros((locations.shape[0], locations.shape[1]))

    def reset(self):
        self.clean()

    def add_value_constraints(self, w: float = 1.0):
        return None

    def add_gradient_constraints(self, w: float = 1.0):
        return None

    def add_norm_constraints(self, w: float = 1.0):
        return None

    def add_tangent_constraints(self, w: float = 1.0):
        return None

    def add_interface_constraints(self, w: float = 1.0):
        return None

    def add_value_inequality_constraints(self, w: float = 1.0):
        return None

    def add_inequality_pairs_constraints(
        self,
        w: float = 1.0,
        upper_bound=np.finfo(float).eps,
        lower_bound=-np.inf,
        pairs=None,
    ):
        return None


def test_default_surfaces_raises_not_implemented():
    interpolator = _MinimalGeologicalInterpolator()

    with pytest.raises(NotImplementedError, match="Surface extraction not implemented"):
        interpolator.surfaces(0.0)


def test_set_constraints_from_pydantic_models(interpolator, data):
    value_rows = data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    normal_rows = data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()

    value_constraint = ValueConstraint.from_array(value_rows)
    normal_constraint = GradientConstraint.from_array(normal_rows, is_normal=True)

    interpolator.set_value_constraints(value_constraint)
    interpolator.set_normal_constraints(normal_constraint)

    assert interpolator.get_value_constraints().shape[0] == value_rows.shape[0]
    assert interpolator.get_norm_constraints().shape[0] == normal_rows.shape[0]


def test_interpolator_json_yaml_round_trip(interpolator, data):
    value_rows = data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    interpolator.set_value_constraints(value_rows)

    json_payload = interpolator.to_json()
    restored_json = GeologicalInterpolator.from_json(json_payload)

    assert restored_json.type == interpolator.type
    assert restored_json.support.n_nodes == interpolator.support.n_nodes
    assert np.array_equal(
        restored_json.get_value_constraints(), interpolator.get_value_constraints()
    )

    yaml_payload = interpolator.to_yaml()
    restored_yaml = GeologicalInterpolator.from_yaml(yaml_payload)

    assert restored_yaml.type == interpolator.type
    assert restored_yaml.support.n_nodes == interpolator.support.n_nodes
    assert np.array_equal(
        restored_yaml.get_value_constraints(), interpolator.get_value_constraints()
    )
