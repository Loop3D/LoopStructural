"""Tests for the opt-in matrix-free StructuredGrid regularisation path.

These tests verify that the ``LinearOperator`` built by
``FiniteDifferenceInterpolator`` when ``regularisation_matrix_free=True`` is
numerically identical (matvec and rmatvec) to the explicit ``coo_matrix`` that
``_assemble_operator``/``add_constraints_to_least_squares`` would otherwise
build, for each of the six second-derivative stencil families (dxx, dyy, dzz,
dxy, dxz, dyz) and the six border first-derivative families
(dx_lower/upper, dy_lower/upper, dz_lower/upper).

The matrix-free path is wired into ``solve_system()``/``fit()`` for both
``cg`` and ``lsmr``. For ``lsmr`` (and any other non-``cg`` solver), the
combined rectangular ``LinearOperator`` from ``get_regularisation_linear_operator``
is used, unchanged. For ``cg``, ``DiscreteInterpolator._solve_with_cg_fused_regularisation``
instead assembles the normal-equations system directly, using
``FiniteDifferenceInterpolator._build_fused_cg_regularisation_operator``'s fused
single-kernel + boundary-corrected regularisation contribution rather than
squaring the combined rectangular operator generically - see the tests below
that exercise this operator directly (``test_fused_cg_regularisation_operator_*``)
as well as the end-to-end solve comparison.
"""

import numpy as np
import pytest
from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid
from loop_interpolation._operator import Operator
from scipy import sparse
from scipy.sparse.linalg import LinearOperator

INTERIOR_OPERATORS = {
    "dxx": Operator.Dxx_mask,
    "dyy": Operator.Dyy_mask,
    "dzz": Operator.Dzz_mask,
    "dxy": Operator.Dxy_mask,
    "dxz": Operator.Dxz_mask,
    "dyz": Operator.Dyz_mask,
}

BORDER_FAMILIES = [
    "dx_lower",
    "dx_upper",
    "dy_lower",
    "dy_upper",
    "dz_lower",
    "dz_upper",
]

ALL_REGULARISATION_FAMILIES = list(INTERIOR_OPERATORS) + BORDER_FAMILIES


def _make_grid() -> StructuredGrid:
    return StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([5, 5, 5]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )


def _explicit_family_matrix(name: str):
    """Assemble a single interior family the normal (explicit coo_matrix) way."""
    interp = FiniteDifferenceInterpolator(_make_grid())
    interp.reset()
    interp._assemble_operator(INTERIOR_OPERATORS[name], 1.0, name=name)
    constraint = interp.constraints[name]
    return constraint["matrix"], constraint["w"], interp.dof


def _matrix_free_family_operator(name: str):
    """Assemble a single interior family via the matrix-free path."""
    interp = FiniteDifferenceInterpolator(_make_grid())
    interp.reset()
    interp.regularisation_matrix_free = True
    interp._assemble_operator(INTERIOR_OPERATORS[name], 1.0, name=name)
    assert name not in interp.constraints, "matrix-free path must not build a coo_matrix"
    op = interp.get_regularisation_linear_operator(names=[name])
    return op, interp.dof


def _explicit_borders():
    interp = FiniteDifferenceInterpolator(_make_grid())
    interp.reset()
    interp.assemble_borders()
    return interp


def _matrix_free_borders():
    interp = FiniteDifferenceInterpolator(_make_grid())
    interp.reset()
    interp.regularisation_matrix_free = True
    interp.assemble_borders()
    return interp


@pytest.mark.parametrize("name", sorted(INTERIOR_OPERATORS))
def test_interior_family_matvec_matches_explicit(name):
    matrix, w, dof = _explicit_family_matrix(name)
    op, mf_dof = _matrix_free_family_operator(name)
    assert mf_dof == dof

    weighted = matrix.multiply(w[:, None]).tocsr()
    rng = np.random.default_rng(0)
    x = rng.normal(size=dof)

    expected = np.asarray(weighted @ x).reshape(-1)
    actual = op.matvec(x)
    diff = np.max(np.abs(actual - expected))
    assert diff < 1e-10, f"{name}: matvec max abs diff {diff}"


@pytest.mark.parametrize("name", sorted(INTERIOR_OPERATORS))
def test_interior_family_rmatvec_matches_explicit(name):
    matrix, w, _dof = _explicit_family_matrix(name)
    op, _ = _matrix_free_family_operator(name)

    weighted = matrix.multiply(w[:, None]).tocsr()
    rng = np.random.default_rng(1)
    y = rng.normal(size=weighted.shape[0])

    expected = np.asarray(weighted.T @ y).reshape(-1)
    actual = op.rmatvec(y)
    diff = np.max(np.abs(actual - expected))
    assert diff < 1e-10, f"{name}: rmatvec max abs diff {diff}"


@pytest.mark.parametrize("name", BORDER_FAMILIES)
def test_border_family_matvec_and_rmatvec_match_explicit(name):
    explicit_interp = _explicit_borders()
    assert name in explicit_interp.constraints
    matrix = explicit_interp.constraints[name]["matrix"]
    w = explicit_interp.constraints[name]["w"]
    weighted = matrix.multiply(w[:, None]).tocsr()

    mf_interp = _matrix_free_borders()
    assert name not in mf_interp.constraints, "matrix-free path must not build a coo_matrix"
    op = mf_interp.get_regularisation_linear_operator(names=[name])
    assert op is not None
    assert op.shape == weighted.shape

    seed = abs(hash(name)) % (2**31)
    rng = np.random.default_rng(seed)
    x = rng.normal(size=mf_interp.dof)
    y = rng.normal(size=weighted.shape[0])

    matvec_diff = np.max(np.abs(op.matvec(x) - np.asarray(weighted @ x).reshape(-1)))
    rmatvec_diff = np.max(np.abs(op.rmatvec(y) - np.asarray(weighted.T @ y).reshape(-1)))
    assert matvec_diff < 1e-10, f"{name}: matvec max abs diff {matvec_diff}"
    assert rmatvec_diff < 1e-10, f"{name}: rmatvec max abs diff {rmatvec_diff}"


def test_combined_operator_matches_full_explicit_regularisation_system():
    """A single combined LinearOperator over all 12 families must match the
    vstack of every family's explicit weighted matrix, for both matvec and
    rmatvec, using non-trivial (non-uniform) per-family weights."""
    setup_kwargs = dict(
        dxx=1.0,
        dyy=0.8,
        dzz=1.3,
        dxy=0.5,
        dxz=0.4,
        dyz=0.6,
        dx=0.7,
        dy=0.9,
        dz=1.1,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )

    explicit_interp = FiniteDifferenceInterpolator(_make_grid())
    explicit_interp.setup_interpolator(**setup_kwargs)

    mats = []
    for name in ALL_REGULARISATION_FAMILIES:
        c = explicit_interp.constraints[name]
        mats.append(c["matrix"].multiply(c["w"][:, None]))
    explicit_reg = sparse.vstack(mats).tocsr()

    mf_interp = FiniteDifferenceInterpolator(_make_grid())
    mf_interp.apply_scaling_matrix = False
    mf_interp.setup_interpolator(regularisation_matrix_free=True, **setup_kwargs)

    for name in ALL_REGULARISATION_FAMILIES:
        assert name not in mf_interp.constraints
        assert name in mf_interp.matrix_free_regularisation_blocks

    op = mf_interp.get_regularisation_linear_operator(names=ALL_REGULARISATION_FAMILIES)
    assert op.shape == explicit_reg.shape

    rng = np.random.default_rng(42)
    x = rng.normal(size=mf_interp.dof)
    y = rng.normal(size=op.shape[0])

    matvec_diff = np.max(np.abs(op.matvec(x) - np.asarray(explicit_reg @ x).reshape(-1)))
    rmatvec_diff = np.max(np.abs(op.rmatvec(y) - np.asarray(explicit_reg.T @ y).reshape(-1)))
    assert matvec_diff < 1e-8, f"combined matvec max abs diff {matvec_diff}"
    assert rmatvec_diff < 1e-8, f"combined rmatvec max abs diff {rmatvec_diff}"


def test_setup_interpolator_matrix_free_flag_skips_explicit_regularisation_matrices():
    interp = FiniteDifferenceInterpolator(_make_grid())
    interp.apply_scaling_matrix = False
    interp.setup_interpolator(
        dxx=1.0,
        dyy=1.0,
        dzz=1.0,
        dxy=1.0,
        dyz=1.0,
        dxz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
        regularisation_matrix_free=True,
    )

    assert interp.regularisation_matrix_free is True
    for name in ALL_REGULARISATION_FAMILIES:
        assert name not in interp.constraints
        assert name in interp.matrix_free_regularisation_blocks

    op = interp.get_regularisation_linear_operator()
    assert op is not None
    assert op.shape[1] == interp.dof


def test_setup_interpolator_matrix_free_flag_falls_back_when_column_scaling_requested():
    interp = FiniteDifferenceInterpolator(_make_grid())
    assert interp.apply_scaling_matrix is True  # default

    interp.setup_interpolator(
        dxx=1.0,
        dyy=1.0,
        dzz=1.0,
        dxy=1.0,
        dyz=1.0,
        dxz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
        regularisation_matrix_free=True,
    )

    # Falls back to the explicit path: flag is reset to False and the usual
    # coo_matrix families are present, nothing recorded as matrix-free.
    assert interp.regularisation_matrix_free is False
    for name in ALL_REGULARISATION_FAMILIES:
        assert name in interp.constraints
    assert interp.matrix_free_regularisation_blocks == {}


def test_regularisation_matrix_free_defaults_to_false():
    interp = FiniteDifferenceInterpolator(_make_grid())
    assert interp.regularisation_matrix_free is False
    assert interp.matrix_free_regularisation_blocks == {}
    interp.setup_interpolator(
        dxx=1.0,
        dyy=1.0,
        dzz=1.0,
        dxy=1.0,
        dyz=1.0,
        dxz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    assert interp.regularisation_matrix_free is False
    for name in ALL_REGULARISATION_FAMILIES:
        assert name in interp.constraints


# ---------------------------------------------------------------------------
# End-to-end wiring: regularisation_matrix_free=True must produce (within
# iterative-solver tolerance) the same solved coefficients / evaluated field
# as the explicit path, for a real solve with value+gradient data constraints
# and nonzero interior regularisation weights.
# ---------------------------------------------------------------------------


def _end_to_end_setup_kwargs():
    return dict(
        dxx=0.4,
        dyy=0.4,
        dzz=0.4,
        dxy=0.1,
        dxz=0.1,
        dyz=0.1,
        dx=0.0,
        dy=0.0,
        dz=0.0,
        cpw=1.0,
        gpw=1.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )


def _build_end_to_end_interpolator(matrix_free: bool) -> FiniteDifferenceInterpolator:
    """Build an FDI with real value + gradient data constraints and nonzero
    dxx/dyy/dzz (+ dxy/dxz/dyz) regularisation, either via the explicit path
    or the matrix-free path (same data, same weights, same solver inputs)."""
    interp = FiniteDifferenceInterpolator(_make_grid())

    rng = np.random.default_rng(7)
    xyz = rng.uniform(1.0, 4.0, size=(12, 3))
    values = np.sin(xyz[:, 0]) + 0.5 * np.cos(xyz[:, 1])
    value_points = np.column_stack([xyz, values, np.ones(xyz.shape[0])])
    interp.set_value_constraints(value_points)

    grad_xyz = rng.uniform(1.0, 4.0, size=(4, 3))
    grad_vectors = np.tile(np.array([0.1, 0.2, 0.97]), (4, 1))
    grad_points = np.column_stack([grad_xyz, grad_vectors, np.ones(4)])
    interp.set_gradient_constraints(grad_points)

    # apply_scaling_matrix=True (the default) forces the matrix-free flag back
    # to the explicit path (see setup_interpolator); use False for both builds
    # so the comparison is solver-choice/regularisation-path only.
    interp.apply_scaling_matrix = False

    interp.setup_interpolator(
        **_end_to_end_setup_kwargs(),
        regularisation_matrix_free=matrix_free,
    )
    return interp


@pytest.mark.parametrize("solver", ["cg", "lsmr"])
def test_matrix_free_regularisation_solve_matches_explicit(solver):
    """The whole point of wiring regularisation_matrix_free into solve_system:
    solving the same problem through the explicit sparse path and through the
    combined LinearOperator path must agree, not just the standalone operator
    matvec/rmatvec checked above."""
    if solver == "cg":
        solver_kwargs = {"maxiter": 5000, "atol": 1e-12, "rtol": 1e-12}
    else:
        solver_kwargs = {"maxiter": 5000, "atol": 1e-12, "btol": 1e-12}

    explicit = _build_end_to_end_interpolator(matrix_free=False)
    matrix_free = _build_end_to_end_interpolator(matrix_free=True)

    assert explicit.regularisation_matrix_free is False
    assert matrix_free.regularisation_matrix_free is True
    assert matrix_free.matrix_free_regularisation_blocks  # non-empty: path actually engaged

    ok_explicit = explicit.solve_system(solver, solver_kwargs=dict(solver_kwargs))
    ok_matrix_free = matrix_free.solve_system(solver, solver_kwargs=dict(solver_kwargs))
    assert ok_explicit is True
    assert ok_matrix_free is True

    # Iterative solvers (cg/lsmr) only converge to within their tolerance, so
    # this is np.allclose (not exact equality). 1e-4 is generous relative to
    # the solved coefficient magnitudes (O(1)) and the requested solver
    # tolerances (1e-12); observed diffs in practice are ~1e-8 (lsmr) to
    # ~1e-11 (cg) for this problem size.
    coeff_diff = np.max(np.abs(explicit.c - matrix_free.c))
    assert coeff_diff < 1e-4, f"{solver}: max abs coefficient diff {coeff_diff}"
    assert np.allclose(explicit.c, matrix_free.c, atol=1e-4, rtol=1e-4)

    sample_points = np.array(
        [
            [2.0, 2.0, 2.0],
            [1.5, 2.5, 1.0],
            [3.0, 1.0, 3.0],
            [2.5, 2.5, 2.5],
        ]
    )
    value_explicit = explicit.evaluate_value(sample_points)
    value_matrix_free = matrix_free.evaluate_value(sample_points)
    value_diff = np.max(np.abs(value_explicit - value_matrix_free))
    assert value_diff < 1e-4, f"{solver}: max abs evaluate_value diff {value_diff}"

    gradient_explicit = explicit.evaluate_gradient(sample_points)
    gradient_matrix_free = matrix_free.evaluate_gradient(sample_points)
    gradient_diff = np.max(np.abs(gradient_explicit - gradient_matrix_free))
    assert gradient_diff < 1e-3, f"{solver}: max abs evaluate_gradient diff {gradient_diff}"


def test_matrix_free_regularisation_falls_back_to_explicit_for_admm_solver():
    """ADMM is out of scope for the matrix-free wiring (it needs an explicit
    sparse system matrix for its inequality-constrained inner solve); selecting
    solver='admm' while regularisation_matrix_free=True must fall back to the
    explicit assembly for that solve rather than crash or silently drop the
    regularisation terms, mirroring the existing apply_scaling_matrix
    fallback pattern.

    Note: this interpolator has no inequality constraints, and ADMM without
    any (a pre-existing, unrelated limitation of the embedded ADMM solver,
    reproducible with regularisation_matrix_free=False too) fails with
    "nelements must be greater than 0" regardless of this feature. So this
    test only asserts the fallback materialisation itself happened - not that
    the ADMM solve succeeded, which is out of scope here.
    """
    interp = _build_end_to_end_interpolator(matrix_free=True)
    assert interp.regularisation_matrix_free is True
    assert interp.matrix_free_regularisation_blocks

    interp.solve_system("admm")

    # solve_system materialised the matrix-free blocks back into explicit
    # constraints before assembling, so nothing matrix-free remains, whether
    # or not the ADMM solve itself went on to succeed.
    assert interp.regularisation_matrix_free is False
    assert interp.matrix_free_regularisation_blocks == {}
    for name in ("dxx", "dyy", "dzz", "dxy", "dxz", "dyz"):
        assert name in interp.constraints


# ---------------------------------------------------------------------------
# Fused single-kernel + boundary-corrected CG regularisation operator
# (`FiniteDifferenceInterpolator._build_fused_cg_regularisation_operator`).
#
# The whole point of this operator is that it must reproduce
# `R_reg^T @ R_reg @ x` EXACTLY (not just in the deep interior, away from the
# true grid boundary) - a fused single (5,5,5)-kernel convolution alone is
# provably wrong in the outer 2-cell shell (composing two radius-1 stencils
# assumes translation invariance, which breaks where the real computation
# discards "rows" that don't correspond to a genuine interior grid node), so
# these tests compare against every single dof, not a masked subset.
# ---------------------------------------------------------------------------


WEIGHTS_DEFAULT = dict(dxx=1.0, dyy=0.8, dzz=1.3, dxy=0.5, dxz=0.4, dyz=0.6)


def _explicit_interior_gram(nsteps, weights=WEIGHTS_DEFAULT):
    """Ground truth: R_reg^T @ R_reg for the 6 interior families only, formed
    as an explicit sparse matrix product - independent of any matrix-free
    machinery."""
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array(nsteps),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )
    interp = FiniteDifferenceInterpolator(grid)
    interp.reset()
    for name, mask in INTERIOR_OPERATORS.items():
        interp._assemble_operator(mask, weights[name], name=name)
    mats = [
        interp.constraints[name]["matrix"].multiply(interp.constraints[name]["w"][:, None])
        for name in INTERIOR_OPERATORS
    ]
    R = sparse.vstack(mats).tocsr()
    return R, interp.dof


def _fused_operator_interior_only(nsteps, weights=WEIGHTS_DEFAULT):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array(nsteps),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )
    interp = FiniteDifferenceInterpolator(grid)
    interp.reset()
    interp.regularisation_matrix_free = True
    for name, mask in INTERIOR_OPERATORS.items():
        interp._assemble_operator(mask, weights[name], name=name)
    op = interp._build_fused_cg_regularisation_operator()
    return op, interp


@pytest.mark.parametrize(
    "nsteps",
    [
        (5, 5, 5),  # matches the other tests in this module; also tiny (max
        # interior distance-to-edge is 2), so almost every dof is in the
        # "boundary shell" this operator must get exactly right.
        (8, 6, 12),  # non-cubic
        (14, 14, 14),
    ],
)
def test_fused_cg_regularisation_operator_matches_explicit_gram_everywhere(nsteps):
    R, dof = _explicit_interior_gram(nsteps)
    op, _interp = _fused_operator_interior_only(nsteps)
    assert op is not None
    assert op.shape == (dof, dof)

    rng = np.random.default_rng(0)
    x = rng.normal(size=dof)
    truth = np.asarray(R.T @ (R @ x)).reshape(-1)
    got = op.matvec(x)

    diff = np.max(np.abs(got - truth))
    assert diff < 1e-10, f"nsteps={nsteps}: max abs diff vs explicit Gram (EVERY dof) {diff}"

    # self-adjoint by construction
    assert np.allclose(got, op.rmatvec(x))


@pytest.mark.parametrize("nsteps", [(6, 6, 6), (9, 7, 5)])
def test_fused_cg_regularisation_operator_is_symmetric(nsteps):
    """<op x, y> == <x, op y> for random x, y: the operator must represent a
    genuinely symmetric matrix (it is a sum of A^T A contributions), not just
    happen to satisfy matvec == rmatvec as a coincidence of implementation."""
    op, interp = _fused_operator_interior_only(nsteps)
    rng = np.random.default_rng(3)
    x = rng.normal(size=interp.dof)
    y = rng.normal(size=interp.dof)
    lhs = np.dot(op.matvec(x), y)
    rhs = np.dot(x, op.matvec(y))
    assert abs(lhs - rhs) < 1e-8 * max(1.0, abs(lhs))


def test_fused_cg_regularisation_operator_matches_explicit_gram_with_border_families():
    """Border first-derivative families (dx_lower/upper, ...) are never part
    of the fused six-family kernel; they must instead be routed through the
    (always-correct, unrestricted) two-pass path unchanged. Verify the full
    12-family combined operator (interior + borders) still matches the
    explicit Gram matrix everywhere."""
    setup_kwargs = dict(
        dxx=1.0,
        dyy=0.8,
        dzz=1.3,
        dxy=0.5,
        dxz=0.4,
        dyz=0.6,
        dx=0.7,
        dy=0.9,
        dz=1.1,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )

    explicit_interp = FiniteDifferenceInterpolator(_make_grid())
    explicit_interp.setup_interpolator(**setup_kwargs)
    mats = []
    for name in ALL_REGULARISATION_FAMILIES:
        c = explicit_interp.constraints[name]
        mats.append(c["matrix"].multiply(c["w"][:, None]))
    R = sparse.vstack(mats).tocsr()

    mf_interp = FiniteDifferenceInterpolator(_make_grid())
    mf_interp.apply_scaling_matrix = False
    mf_interp.setup_interpolator(regularisation_matrix_free=True, **setup_kwargs)

    op = mf_interp._build_fused_cg_regularisation_operator()
    assert op is not None
    assert op.shape == (mf_interp.dof, mf_interp.dof)

    rng = np.random.default_rng(11)
    x = rng.normal(size=mf_interp.dof)
    truth = np.asarray(R.T @ (R @ x)).reshape(-1)
    got = op.matvec(x)
    diff = np.max(np.abs(got - truth))
    assert diff < 1e-9, f"max abs diff vs explicit 12-family Gram (EVERY dof) {diff}"


def test_fused_cg_regularisation_operator_falls_back_safely_for_nonuniform_weights():
    """`use_regularisation_weight_scale=True` makes per-row regularisation
    weight vary spatially, which breaks the translation-invariance assumption
    the fused kernel relies on. The six interior families must then be routed
    to the exact (unfused) two-pass path instead of being fused - never
    silently wrong. Verify this still matches the explicit Gram everywhere."""
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([9, 7, 11]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )

    def _build(matrix_free):
        interp = FiniteDifferenceInterpolator(grid)
        interp.reset()
        interp.use_regularisation_weight_scale = True
        interp.regularisation_scale = np.linspace(0.5, 2.0, interp.dof)
        interp.regularisation_matrix_free = matrix_free
        for name, mask in INTERIOR_OPERATORS.items():
            interp._assemble_operator(mask, WEIGHTS_DEFAULT[name], name=name)
        return interp

    explicit_interp = _build(False)
    mats = [
        explicit_interp.constraints[name]["matrix"].multiply(
            explicit_interp.constraints[name]["w"][:, None]
        )
        for name in INTERIOR_OPERATORS
    ]
    R = sparse.vstack(mats).tocsr()

    mf_interp = _build(True)
    for name in INTERIOR_OPERATORS:
        assert name in mf_interp.matrix_free_regularisation_blocks
        w = mf_interp.matrix_free_regularisation_blocks[name]["w"]
        assert not np.all(w == w[0]), "test setup should produce non-uniform weights"

    op = mf_interp._build_fused_cg_regularisation_operator()
    assert op is not None

    rng = np.random.default_rng(5)
    x = rng.normal(size=mf_interp.dof)
    truth = np.asarray(R.T @ (R @ x)).reshape(-1)
    got = op.matvec(x)
    diff = np.max(np.abs(got - truth))
    assert diff < 1e-9, f"max abs diff vs explicit Gram with non-uniform weights {diff}"


def test_use_fused_cg_regularisation_predicate():
    interp = _build_end_to_end_interpolator(matrix_free=True)
    assert interp._use_fused_cg_regularisation("cg") is True
    assert interp._use_fused_cg_regularisation("lsmr") is False
    assert interp._use_fused_cg_regularisation("admm") is False

    explicit_interp = _build_end_to_end_interpolator(matrix_free=False)
    assert explicit_interp._use_fused_cg_regularisation("cg") is False


def test_solve_system_cg_actually_takes_fused_path():
    """Sanity that the end-to-end cg solve test above is really exercising the
    new fused-hybrid internals, not silently falling back."""
    interp = _build_end_to_end_interpolator(matrix_free=True)
    assert interp._use_fused_cg_regularisation("cg") is True
    ok = interp.solve_system("cg", solver_kwargs={"maxiter": 5000, "atol": 1e-12, "rtol": 1e-12})
    assert ok is True
    # matrix-free state must still be intact afterwards (cg fast path doesn't
    # materialise/clear it the way the admm fallback does)
    assert interp.regularisation_matrix_free is True
    assert interp.matrix_free_regularisation_blocks


def test_lsmr_still_uses_combined_rectangular_linear_operator():
    """lsmr must be completely unaffected by the cg fast path: build_matrix()
    still returns the combined rectangular LinearOperator (data rows stacked
    on regularisation rows) for it, exactly as before."""
    interp = _build_end_to_end_interpolator(matrix_free=True)
    assert interp._use_fused_cg_regularisation("lsmr") is False

    A, b = interp.build_matrix()
    assert isinstance(A, LinearOperator)
    n_data = sum(len(c["w"]) for c in interp.constraints.values())
    n_reg = sum(block["idc"].shape[0] for block in interp.matrix_free_regularisation_blocks.values())
    assert A.shape == (n_data + n_reg, interp.dof)
    assert b.shape[0] == n_data + n_reg

    ok = interp.solve_system("lsmr", solver_kwargs={"maxiter": 5000, "atol": 1e-12, "btol": 1e-12})
    assert ok is True
