import numpy as np
import pytest
from loop_interpolation import (
    DirectionalRegularisation,
    FiniteDifferenceInterpolator,
    PiecewiseLinearInterpolator,
    RegularisationConfig,
    StructuredGrid,
    TetMesh,
    add_element_anisotropy_constraints,
)


def _constant_direction(points: np.ndarray, direction=(0.0, 0.0, 1.0)) -> np.ndarray:
    vector = np.asarray(direction, dtype=float)
    vector /= np.linalg.norm(vector)
    return np.tile(vector, (points.shape[0], 1))


def _make_structured_grid() -> StructuredGrid:
    return StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([5, 5, 5]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )


def _make_tet_mesh() -> TetMesh:
    return TetMesh(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([5, 5, 5]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )


@pytest.mark.parametrize(
    ("factory", "setup_kwargs"),
    (
        (
            lambda: FiniteDifferenceInterpolator(_make_structured_grid()),
            {"dxx": 0.0, "dyy": 0.0, "dzz": 0.0, "dxy": 0.0, "dyz": 0.0, "dxz": 0.0},
        ),
        (lambda: PiecewiseLinearInterpolator(_make_tet_mesh()), {"cgw": 0.0}),
    ),
)
def test_shared_directional_regularisation_dict_works_across_support_types(factory, setup_kwargs):
    interpolator = factory()
    regularisation = {
        "isotropic": 0.0,
        "directional": [
            {
                "weight": 2.5,
                "direction": lambda points: _constant_direction(points, direction=(0.0, 0.0, 1.0)),
                "name": "shared vertical smoothing",
            }
        ],
    }

    interpolator.setup_interpolator(
        regularisation=regularisation,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
        **setup_kwargs,
    )

    matching = [
        name for name in interpolator.constraints if name.startswith("shared vertical smoothing")
    ]
    assert matching
    assert all(interpolator.constraints[name]["matrix"].shape[0] > 0 for name in matching)


def test_shared_directional_regularisation_config_object_is_accepted():
    interpolator = PiecewiseLinearInterpolator(_make_tet_mesh())
    regularisation = RegularisationConfig(
        isotropic=0.0,
        directional=(
            DirectionalRegularisation(
                weight=1.0,
                direction=lambda points: _constant_direction(points, direction=(1.0, 0.0, 0.0)),
                name="config object smoothing",
            ),
        ),
    )

    interpolator.setup_interpolator(
        regularisation=regularisation,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
        cgw=0.0,
    )

    assert any(name.startswith("config object smoothing") for name in interpolator.constraints)


def test_anisotropic_p1_regularisation_uses_shared_directional_api():
    class _AnisotropyStub:
        def get_directions(self, points):
            primary = _constant_direction(points, direction=(1.0, 0.0, 0.0))
            secondary = _constant_direction(points, direction=(0.0, 1.0, 0.0))
            normal = _constant_direction(points, direction=(0.0, 0.0, 1.0))
            return primary, secondary, normal

    interpolator = PiecewiseLinearInterpolator(_make_tet_mesh())
    interpolator.setup_interpolator(
        cgw=0.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    add_element_anisotropy_constraints(
        interpolator,
        _AnisotropyStub(),
        primary_w=None,
        secondary_w=None,
        normal_w=None,
        regularisation_weights=[0.1, 0.01, 0.01],
    )

    matching = [name for name in interpolator.constraints if "anisotropy regularisation" in name]
    assert matching
    assert all(interpolator.constraints[name]["matrix"].shape[0] > 0 for name in matching)


def test_p1_regularisation_weight_scale_creates_spatially_varying_weights():
    interpolator = PiecewiseLinearInterpolator(_make_tet_mesh())
    normal_constraints = np.array([[2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 1.0]])
    interpolator.set_normal_constraints(normal_constraints)

    interpolator.setup_interpolator(
        cgw=0.1,
        cpw=0.0,
        gpw=0.0,
        npw=1.0,
        tpw=0.0,
        ipw=0.0,
        use_regularisation_weight_scale=True,
    )

    edge_jump = interpolator.constraints["edge jump"]
    assert edge_jump["w"].size > 1
    assert np.ptp(edge_jump["w"]) > 0.0


def test_p1_regularisation_weight_sigma_controls_decay_strength():
    normal_constraints = np.array([[2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 1.0]])

    local = PiecewiseLinearInterpolator(_make_tet_mesh())
    local.set_normal_constraints(normal_constraints)
    local.setup_interpolator(
        cgw=0.1,
        cpw=0.0,
        gpw=0.0,
        npw=1.0,
        tpw=0.0,
        ipw=0.0,
        use_regularisation_weight_scale=True,
        regularisation_weight_sigma=0.2,
    )

    broad = PiecewiseLinearInterpolator(_make_tet_mesh())
    broad.set_normal_constraints(normal_constraints)
    broad.setup_interpolator(
        cgw=0.1,
        cpw=0.0,
        gpw=0.0,
        npw=1.0,
        tpw=0.0,
        ipw=0.0,
        use_regularisation_weight_scale=True,
        regularisation_weight_sigma=2.0,
    )

    assert np.mean(broad.constraints["edge jump"]["w"]) > np.mean(
        local.constraints["edge jump"]["w"]
    )


def test_fdi_regularisation_weight_sigma_controls_decay_strength():
    normal_constraints = np.array([[2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 1.0]])

    local = FiniteDifferenceInterpolator(_make_structured_grid())
    local.set_normal_constraints(normal_constraints)
    local.setup_interpolator(
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        cpw=0.0,
        gpw=0.0,
        npw=1.0,
        tpw=0.0,
        ipw=0.0,
        use_regularisation_weight_scale=True,
        regularisation_weight_sigma=0.2,
    )

    broad = FiniteDifferenceInterpolator(_make_structured_grid())
    broad.set_normal_constraints(normal_constraints)
    broad.setup_interpolator(
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        cpw=0.0,
        gpw=0.0,
        npw=1.0,
        tpw=0.0,
        ipw=0.0,
        use_regularisation_weight_scale=True,
        regularisation_weight_sigma=2.0,
    )

    assert np.ptp(local.regularisation_scale) > np.ptp(broad.regularisation_scale)
