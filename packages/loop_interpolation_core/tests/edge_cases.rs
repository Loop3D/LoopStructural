//! Adversarial / edge-case coverage for `loop_interpolation_core`.
//!
//! Unlike the inline `#[cfg(test)]` modules in `src/*.rs` (which exercise
//! the happy path: normal-operator symmetry, affine-field recovery,
//! stencil/CG correctness), this file targets degenerate or malformed
//! inputs: single-node grid axes, zero/negative regularisation weights,
//! NaN/Inf constraint and query data, and extreme step-vector anisotropy.
//! The goal throughout is "does this panic / silently corrupt state", not
//! "is this numerically ideal" -- several of these document existing,
//! acceptable behaviour (e.g. NaN propagating into a solve) rather than
//! asserting an idealised response.
//!
//! This is a Cargo *integration* test: it only uses the crate's public API
//! (`loop_interpolation_core::{...}`), deliberately kept disjoint from the
//! new constraint-type work landing concurrently on another
//! branch in `src/grid.rs` / `src/stencil.rs` / `src/interpolator.rs`.

use loop_interpolation_core::interpolator::BorderRegularisationWeights;
use loop_interpolation_core::stencil::{dxx, dxy, dxz, dyy, dyz, dzz};
use loop_interpolation_core::{
    FiniteDifferenceInterpolator, RegularisationWeights, StructuredGrid3D,
};

fn all_finite(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite())
}

// ---------------------------------------------------------------------
// Degenerate grids: a `1` in one or more `nsteps` axes.
// ---------------------------------------------------------------------

#[test]
fn trilinear_coeffs_none_for_single_node_axis_x() {
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [1, 5, 5], [1.0, 1.0, 1.0]);
    // Even a point that is otherwise perfectly "inside" the y/z extent
    // must come back None because the x axis has < 2 nodes.
    assert!(grid.trilinear_coeffs([0.0, 2.0, 2.0]).is_none());
    assert!(grid.trilinear_coeffs([0.0, 0.0, 0.0]).is_none());
}

#[test]
fn trilinear_coeffs_none_for_single_node_axis_y_and_z() {
    let grid_y = StructuredGrid3D::new([0.0, 0.0, 0.0], [5, 1, 5], [1.0, 1.0, 1.0]);
    assert!(grid_y.trilinear_coeffs([2.0, 0.0, 2.0]).is_none());

    let grid_z = StructuredGrid3D::new([0.0, 0.0, 0.0], [5, 5, 1], [1.0, 1.0, 1.0]);
    assert!(grid_z.trilinear_coeffs([2.0, 2.0, 0.0]).is_none());
}

#[test]
fn trilinear_coeffs_none_for_fully_degenerate_grid() {
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [1, 1, 1], [1.0, 1.0, 1.0]);
    assert_eq!(grid.n_nodes(), 1);
    assert!(grid.trilinear_coeffs([0.0, 0.0, 0.0]).is_none());
}

#[test]
fn interpolator_on_single_node_axis_grid_does_not_panic() {
    // nsteps = [1, 5, 5]: every point is rejected by trilinear_coeffs (the
    // degenerate axis short-circuits before the bounding-box check), so
    // constraints are silently dropped -- same contract as out-of-bounds
    // points -- and solve()/evaluate() must still run cleanly on the
    // resulting (constraint-free) system.
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [1, 5, 5], [1.0, 1.0, 1.0]);
    let mut interp = FiniteDifferenceInterpolator::new(grid);
    interp.add_value_constraint([0.0, 2.0, 2.0], 1.0, 5.0);
    assert_eq!(interp.n_value_constraints(), 0);

    assert_eq!(interp.dof(), 5 * 5);
    let diag = interp.diagonal();
    assert_eq!(diag.len(), interp.dof());
    assert!(all_finite(&diag));

    let x = vec![0.0; interp.dof()];
    let n_x = interp.apply_normal_operator(&x);
    assert!(all_finite(&n_x));

    let result = interp.solve(50, 1e-10);
    assert!(
        result.converged,
        "expected trivial all-zero-rhs system to converge immediately"
    );
    assert!(all_finite(&result.x));

    let evaluated = interp.evaluate(&result.x, &[[0.0, 2.0, 2.0]]);
    // Query point is rejected for the same reason value constraints are.
    assert!(evaluated[0].is_nan());
}

#[test]
fn interpolator_on_fully_degenerate_grid_does_not_panic() {
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [1, 1, 1], [1.0, 1.0, 1.0]);
    let mut interp = FiniteDifferenceInterpolator::new(grid);
    interp.add_value_constraint([0.0, 0.0, 0.0], 3.0, 1.0);
    assert_eq!(interp.n_value_constraints(), 0);
    assert_eq!(interp.dof(), 1);

    let result = interp.solve(10, 1e-10);
    assert!(result.converged);
    assert_eq!(result.x.len(), 1);
    assert!(result.x[0].is_finite());
}

// ---------------------------------------------------------------------
// Zero / negative regularisation weights.
// ---------------------------------------------------------------------

fn make_grid_9() -> StructuredGrid3D {
    StructuredGrid3D::new([0.0, 0.0, 0.0], [9, 9, 9], [1.0, 1.0, 1.0])
}

fn sample_field(n: usize, seed: u64) -> Vec<f64> {
    (0..n)
        .map(|i| (((i as u64).wrapping_mul(seed).wrapping_add(11)) % 97) as f64 / 97.0)
        .collect()
}

#[test]
fn zero_regularisation_weights_drop_stencil_contribution() {
    // With every family weight at 0.0, only the data term + ridge damping
    // should contribute to the normal operator/diagonal.
    let grid = make_grid_9();
    let zero_weights = RegularisationWeights {
        dxx: 0.0,
        dyy: 0.0,
        dzz: 0.0,
        dxy: 0.0,
        dxz: 0.0,
        dyz: 0.0,
    };
    // Border (first-derivative) regularisation was added by a later P0 pass
    // and defaults to weight 1.0, always on -- zero it out too, since this
    // test is specifically about the six curvature families contributing
    // nothing, not about border regularisation.
    let mut interp = FiniteDifferenceInterpolator::new(grid)
        .with_regularisation(zero_weights)
        .with_border_regularisation(BorderRegularisationWeights {
            dx: 0.0,
            dy: 0.0,
            dz: 0.0,
        });
    interp.add_value_constraint([4.0, 4.0, 4.0], 1.0, 1.0);

    let x = sample_field(interp.dof(), 31);
    let n_x = interp.apply_normal_operator(&x);
    assert!(all_finite(&n_x));

    // Reference: data-only contribution + ridge, computed by hand from the
    // same public API (ridge_factor defaults to 1e-8).
    let ridge2 = interp.ridge_factor * interp.ridge_factor;
    let expected: Vec<f64> = x.iter().map(|xi| ridge2 * xi).collect();
    // n_x should equal ridge term plus the (small, localised) data-term
    // contribution -- rather than re-derive that by hand, just check that
    // away from the single constrained node's 8-corner footprint, n_x is
    // exactly the ridge term (i.e. no stray curvature contribution leaks
    // in when all six weights are zero).
    let grid = interp.grid;
    let touched: std::collections::HashSet<usize> = grid
        .trilinear_coeffs([4.0, 4.0, 4.0])
        .unwrap()
        .0
        .iter()
        .copied()
        .collect();
    for (idx, (&nxi, &exp)) in n_x.iter().zip(expected.iter()).enumerate() {
        if !touched.contains(&idx) {
            assert!((nxi - exp).abs() < 1e-12, "node {idx}: {nxi} vs {exp}");
        }
    }
}

#[test]
fn negative_regularisation_weights_match_positive_counterpart() {
    // Weights are always squared before use, so -w must behave identically
    // to +w in every observable: apply_normal_operator, diagonal, and the
    // solved field.
    let grid = make_grid_9();
    let positive = RegularisationWeights {
        dxx: 2.0,
        dyy: 1.5,
        dzz: 0.7,
        dxy: 1.0,
        dxz: 1.0,
        dyz: 1.0,
    };
    let negative = RegularisationWeights {
        dxx: -2.0,
        dyy: -1.5,
        dzz: -0.7,
        dxy: -1.0,
        dxz: -1.0,
        dyz: -1.0,
    };

    let build = |weights: RegularisationWeights| {
        let mut interp =
            FiniteDifferenceInterpolator::new(make_grid_9()).with_regularisation(weights);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 3.0);
        interp.add_value_constraint([6.0, 3.0, 4.0], -0.5, 2.0);
        interp
    };
    let _ = grid; // grid var only used to anchor dof() below

    let interp_pos = build(positive);
    let interp_neg = build(negative);

    let x = sample_field(interp_pos.dof(), 17);
    let n_pos = interp_pos.apply_normal_operator(&x);
    let n_neg = interp_neg.apply_normal_operator(&x);
    for (a, b) in n_pos.iter().zip(n_neg.iter()) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }

    let diag_pos = interp_pos.diagonal();
    let diag_neg = interp_neg.diagonal();
    for (a, b) in diag_pos.iter().zip(diag_neg.iter()) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }

    let solved_pos = interp_pos.solve(500, 1e-10);
    let solved_neg = interp_neg.solve(500, 1e-10);
    assert!(solved_pos.converged && solved_neg.converged);
    for (a, b) in solved_pos.x.iter().zip(solved_neg.x.iter()) {
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
    }
}

#[test]
fn negative_ridge_factor_matches_positive_counterpart() {
    // ridge_factor is squared the same way regularisation weights are.
    let mut pos = FiniteDifferenceInterpolator::new(make_grid_9());
    pos.ridge_factor = 0.05;
    pos.add_value_constraint([4.0, 4.0, 4.0], 2.0, 1.0);

    let mut neg = FiniteDifferenceInterpolator::new(make_grid_9());
    neg.ridge_factor = -0.05;
    neg.add_value_constraint([4.0, 4.0, 4.0], 2.0, 1.0);

    let x = sample_field(pos.dof(), 5);
    let n_pos = pos.apply_normal_operator(&x);
    let n_neg = neg.apply_normal_operator(&x);
    for (a, b) in n_pos.iter().zip(n_neg.iter()) {
        assert!((a - b).abs() < 1e-9);
    }
}

#[test]
fn zero_weight_value_constraint_contributes_nothing() {
    let grid = make_grid_9();
    let mut with_extra_zero_weight = FiniteDifferenceInterpolator::new(grid);
    with_extra_zero_weight.add_value_constraint([2.0, 2.0, 2.0], 1.0, 4.0);
    with_extra_zero_weight.add_value_constraint([6.0, 6.0, 6.0], 999.0, 0.0); // zero-weight: should be a no-op

    let baseline_grid = make_grid_9();
    let mut baseline = FiniteDifferenceInterpolator::new(baseline_grid);
    baseline.add_value_constraint([2.0, 2.0, 2.0], 1.0, 4.0);

    // The zero-weight constraint is still recorded...
    assert_eq!(with_extra_zero_weight.n_value_constraints(), 2);
    assert_eq!(baseline.n_value_constraints(), 1);

    // ...but has no effect on diagonal(), rhs(), or the normal operator.
    let diag_a = with_extra_zero_weight.diagonal();
    let diag_b = baseline.diagonal();
    for (a, b) in diag_a.iter().zip(diag_b.iter()) {
        assert!((a - b).abs() < 1e-15, "{a} vs {b}");
    }

    let rhs_a = with_extra_zero_weight.rhs();
    let rhs_b = baseline.rhs();
    for (a, b) in rhs_a.iter().zip(rhs_b.iter()) {
        assert!((a - b).abs() < 1e-15, "{a} vs {b}");
    }

    let x = sample_field(with_extra_zero_weight.dof(), 41);
    let n_a = with_extra_zero_weight.apply_normal_operator(&x);
    let n_b = baseline.apply_normal_operator(&x);
    for (a, b) in n_a.iter().zip(n_b.iter()) {
        assert!((a - b).abs() < 1e-15, "{a} vs {b}");
    }
}

// ---------------------------------------------------------------------
// NaN / Inf inputs.
// ---------------------------------------------------------------------

#[test]
fn nan_value_propagates_into_rhs_but_not_diagonal() {
    // Documents current behaviour: a NaN *value* is not rejected the way
    // an out-of-bounds *point* is -- it flows straight into rhs() (and
    // from there into solve()'s residual), while diagonal() (which never
    // reads `.value`) stays finite. This is a known gap, not a guarantee.
    let mut interp = FiniteDifferenceInterpolator::new(make_grid_9());
    interp.add_value_constraint([4.0, 4.0, 4.0], f64::NAN, 1.0);
    assert_eq!(interp.n_value_constraints(), 1);

    let diag = interp.diagonal();
    assert!(
        all_finite(&diag),
        "diagonal() must not be poisoned by a NaN constraint value"
    );

    let rhs = interp.rhs();
    assert!(
        rhs.iter().any(|v| v.is_nan()),
        "NaN constraint value should surface in rhs()"
    );

    // NaN in one rhs entry poisons the global dot-products CG relies on,
    // so the solved field ends up NaN too -- this is the actual, currently
    // unguarded behaviour, asserted here rather than "fixed" silently.
    let result = interp.solve(20, 1e-10);
    assert!(result.x.iter().any(|v| v.is_nan()));
}

#[test]
fn infinite_weight_produces_infinite_diagonal_entry() {
    // Documents current behaviour: weight is squared and used directly, so
    // an infinite weight yields an infinite diagonal entry rather than
    // being rejected the way a zero-norm row already is.
    let mut interp = FiniteDifferenceInterpolator::new(make_grid_9());
    interp.add_value_constraint([4.0, 4.0, 4.0], 1.0, f64::INFINITY);
    assert_eq!(interp.n_value_constraints(), 1);

    let diag = interp.diagonal();
    assert!(diag.iter().any(|v| v.is_infinite()));
}

#[test]
fn nan_query_point_is_now_dropped_like_out_of_bounds() {
    // Before the fix in `StructuredGrid3D::trilinear_coeffs`, a NaN
    // coordinate slipped past the bounds check (every comparison against
    // NaN is false) and was then silently clamped to 0.0 by `f64::clamp`,
    // i.e. treated as if the point sat exactly at the grid's lower corner
    // -- a real point silently substituted for a garbage one. The bounds
    // check now uses the "reject unless provably inside" form, which
    // recovers `None` for NaN the same way it already did for +-inf.
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [9, 9, 9], [1.0, 1.0, 1.0]);
    assert!(grid.trilinear_coeffs([f64::NAN, 2.0, 2.0]).is_none());
    assert!(grid.trilinear_coeffs([2.0, f64::NAN, 2.0]).is_none());
    assert!(grid.trilinear_coeffs([2.0, 2.0, f64::NAN]).is_none());

    // +-inf were already rejected pre-fix; assert they still are.
    assert!(grid.trilinear_coeffs([f64::INFINITY, 2.0, 2.0]).is_none());
    assert!(
        grid.trilinear_coeffs([f64::NEG_INFINITY, 2.0, 2.0])
            .is_none()
    );

    let mut interp = FiniteDifferenceInterpolator::new(grid);
    interp.add_value_constraint([f64::NAN, 2.0, 2.0], 1.0, 1.0);
    assert_eq!(
        interp.n_value_constraints(),
        0,
        "NaN-coordinate constraint should be dropped, not silently relocated"
    );

    let result = interp.solve(10, 1e-10);
    let evaluated = interp.evaluate(
        &result.x,
        &[[f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]],
    );
    assert!(evaluated[0].is_nan());
    assert!(evaluated[1].is_nan());
}

// ---------------------------------------------------------------------
// Extremely anisotropic step vectors.
// ---------------------------------------------------------------------

#[test]
fn trilinear_coeffs_finite_and_sum_to_one_under_extreme_anisotropy() {
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [5, 5, 5], [1e-6, 1.0, 1e6]);
    // Interior, non-node-aligned point well within each axis's extent
    // (x extent = 4e-6, y extent = 4.0, z extent = 4e6).
    let point = [1.5e-6, 2.5, 1.5e6];
    let (idx, coefs) = grid
        .trilinear_coeffs(point)
        .expect("point is within bounds");
    assert!(all_finite(&coefs));
    let sum: f64 = coefs.iter().sum();
    assert!((sum - 1.0).abs() < 1e-9, "sum={sum}");
    assert!(idx.iter().all(|&i| i < grid.n_nodes()));
}

#[test]
fn curvature_stencils_stay_finite_under_extreme_anisotropy() {
    // The stencils operate purely on grid-index space (they don't read
    // step_vector at all), so this mostly documents that fact -- an
    // extreme step_vector cannot itself introduce overflow/underflow into
    // forward/adjoint, since no division by step size ever happens there.
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [7, 7, 7], [1e-6, 1.0, 1e6]);
    let x: Vec<f64> = (0..grid.n_nodes())
        .map(|flat| {
            let (i, j, k) = grid.node_ijk(flat);
            (i * i + j + 3 * k) as f64
        })
        .collect();

    for stencil in [dxx(), dyy(), dzz(), dxy(), dxz(), dyz()] {
        let y = stencil.forward(&grid, &x);
        assert!(all_finite(&y), "stencil forward produced non-finite output");
        let adj = stencil.adjoint(&grid, &y);
        assert!(
            all_finite(&adj),
            "stencil adjoint produced non-finite output"
        );
    }
}

#[test]
fn solve_stays_finite_under_extreme_anisotropy() {
    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [5, 5, 5], [1e-6, 1.0, 1e6]);
    let mut interp = FiniteDifferenceInterpolator::new(grid);

    // Scatter a handful of value constraints across the (physically tiny
    // in x, huge in z) domain.
    let samples = [
        ([0.5e-6, 0.5, 0.5e6], 1.0),
        ([3.5e-6, 3.5, 3.5e6], 2.0),
        ([1.0e-6, 3.0, 2.0e6], -1.0),
        ([2.0e-6, 1.0, 0.5e6], 0.5),
    ];
    for (p, v) in samples {
        interp.add_value_constraint(p, v, 10.0);
    }
    assert_eq!(interp.n_value_constraints(), samples.len());

    let diag = interp.diagonal();
    assert!(
        all_finite(&diag),
        "diagonal() must stay finite under extreme step-vector anisotropy"
    );

    let result = interp.solve(2000, 1e-8);
    assert!(
        all_finite(&result.x),
        "solve() produced non-finite node values under extreme anisotropy"
    );

    let query = [[2.0e-6, 2.0, 2.0e6]];
    let predicted = interp.evaluate(&result.x, &query);
    assert!(all_finite(&predicted));
}
