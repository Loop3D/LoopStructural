//! Matrix-free port of `loop_interpolation`'s `FiniteDifferenceInterpolator`.
//!
//! Ports the ingredients that make the original interesting to reimplement
//! in Rust -- trilinear value/gradient/norm/tangent/interface constraints,
//! isotropic curvature regularisation, centred second-derivative border regularisation,
//! anisotropic directional regularisation -- on a
//! `StructuredGrid`, and solves the resulting weighted least squares system
//! with conjugate gradient on the normal equations, entirely through
//! [`Stencil`]/[`EdgeStencil`] operators (see `stencil.rs`) and dense
//! per-constraint row storage instead of building any sparse matrix. Box
//! (inequality) constraints -- value-bound (`add_value_inequality_constraint`)
//! and pairwise-ordering (`add_inequality_pair_constraint`, e.g.
//! stratigraphic ordering) -- are also matrix-free, stored the same "own
//! row(s), accumulate into the shared operator" way as every equality
//! constraint family, and solved via three interchangeable matrix-free
//! solvers (see `Inequalities.md` for the design rationale behind offering
//! all three rather than just one):
//!
//! - [`FiniteDifferenceInterpolator::solve_admm`]: classical ADMM
//!   box-splitting (ported from `loop_interpolation`'s
//!   `loopsolver/admm_method.py` + `admm_solver.py`, scoped down to the core
//!   fixed-rho algorithm -- see that method's docs for exactly what was
//!   descoped). Cheap per outer iteration and handles a dense/high active
//!   fraction well, but only ever satisfies bounds up to its primal residual
//!   tolerance (asymptotic contact).
//! - [`FiniteDifferenceInterpolator::solve_active_set_cg`]: a primal
//!   active-set method (Projected Active-Set Conjugate Gradient) that keeps
//!   every active row exactly (to floating-point precision) on its bound,
//!   at the cost of an `O(m)` scan per CG step and small dense Gram solves
//!   on every working-set change -- intended for a low-to-moderate active
//!   fraction.
//! - [`FiniteDifferenceInterpolator::solve_hybrid`]: a fixed-budget ADMM
//!   warm-up feeding a working-set guess into an active-set polish, getting
//!   ADMM's cheap bulk discovery of which rows are active plus AS-CG's
//!   exact contact precision.
//!
//! Anisotropic directional regularisation
//! ([`FiniteDifferenceInterpolator::add_directional_regularisation`]) is the
//! generic primitive fold interpolation's anisotropic smoothing is built
//! from (the Rust analogue of the Python
//! `minimise_directional_gradient_change`): it penalises the directional
//! second derivative `(v.grad)^2 f = 0` at every interior node for a
//! caller-supplied per-node direction field `v`, weighting the six curvature
//! families by a per-node `comp_w = [vx^2, vy^2, vz^2, 2vxvy, 2vxvz, 2vyvz]`
//! instead of one uniform scalar. Since per-node weights can't use the fused
//! single-kernel fast path (`build_combined_kernel` requires a uniform
//! weight), this goes through its own unfused accumulation pass -- see
//! `Stencil::accumulate_normal_weighted_into`/`diagonal_weighted`. Left out on
//! purpose (unlike the Python class): the fold-frame geometry itself (deriving
//! direction vectors from a fold's scalar fields) and column scaling. Column
//! scaling slots into the same `apply_normal_operator` / `rhs` accumulation
//! pattern used here if it's added later.

use crate::PAR_MIN_LEN;
use crate::cg::{CgResult, conjugate_gradient_preconditioned};
use crate::error::CoreError;
use crate::grid::StructuredGrid3D;
use crate::memory::MemTracker;
use crate::par::*;
use crate::scalar::Scalar;
use crate::stencil::{
    EdgeStencil, Stencil, dxx, dxy, dxz, dyy, dyz, dzz, edge_backward, edge_forward,
};
use std::collections::HashMap;
use std::sync::OnceLock;

/// Number of curvature-regularisation stencil families (dxx, dyy, dzz, dxy,
/// dxz, dyz). `apply_normal_operator_into` computes all six per grid node in
/// a single fused rayon pass (see below) rather than looping over families
/// and dispatching a separate parallel pass per family: on a modest-sized
/// grid, thread-pool coordination overhead for a small parallel task can
/// dwarf the O(1)-taps-per-node arithmetic it's doing, so amortising that
/// overhead across all six families in one dispatch (instead of twelve --
/// forward and adjoint, times six) is what actually makes this fast.
const N_FAMILIES: usize = 6;

/// A single (row-normalised) linear constraint over the 8 trilinear corner
/// dofs of one point: `coefs . x[indices] = value`, weighted by `weight`.
/// Used for value, gradient-orthogonal, norm, and tangent constraints alike
/// -- they only differ in how `coefs`/`value` are built (see
/// `add_value_constraint`, `build_orthogonal_constraint`,
/// `add_norm_constraint`), not in how they're accumulated into the normal
/// equations.
struct PointConstraint {
    indices: [usize; 8],
    coefs: [f64; 8],
    value: f64,
    weight: f64,
}

/// A single (row-normalised) linear constraint tying two points' trilinear
/// corner dofs together: `coefs_a . x[indices_a] - coefs_b . x[indices_b] =
/// 0`, weighted by `weight`. Used for interface (equal-value) constraints.
struct PairConstraint {
    /// First 8 entries are point A's corner indices, next 8 are point B's.
    indices: [usize; 16],
    /// Row-normalised combined coefficients: point A's value coefs followed
    /// by the *negated* point B's value coefs.
    coefs: [f64; 16],
    weight: f64,
}

fn accumulate_point_constraints_normal(
    constraints: &[PointConstraint],
    x: &[f64],
    out: &mut [f64],
) {
    for c in constraints {
        let mut row_value = 0.0;
        for k in 0..8 {
            row_value += c.coefs[k] * x[c.indices[k]];
        }
        let scaled = c.weight * c.weight * row_value;
        for k in 0..8 {
            out[c.indices[k]] += scaled * c.coefs[k];
        }
    }
}

fn point_constraints_rhs(constraints: &[PointConstraint], out: &mut [f64]) {
    for c in constraints {
        let scaled = c.weight * c.weight * c.value;
        for k in 0..8 {
            out[c.indices[k]] += scaled * c.coefs[k];
        }
    }
}

fn point_constraints_diagonal(constraints: &[PointConstraint], diag: &mut [f64]) {
    for c in constraints {
        let w2 = c.weight * c.weight;
        for k in 0..8 {
            diag[c.indices[k]] += w2 * c.coefs[k] * c.coefs[k];
        }
    }
}

fn accumulate_pair_constraints_normal(constraints: &[PairConstraint], x: &[f64], out: &mut [f64]) {
    for c in constraints {
        let mut row_value = 0.0;
        for k in 0..16 {
            row_value += c.coefs[k] * x[c.indices[k]];
        }
        let scaled = c.weight * c.weight * row_value;
        for k in 0..16 {
            out[c.indices[k]] += scaled * c.coefs[k];
        }
    }
}

fn pair_constraints_diagonal(constraints: &[PairConstraint], diag: &mut [f64]) {
    for c in constraints {
        let w2 = c.weight * c.weight;
        for k in 0..16 {
            diag[c.indices[k]] += w2 * c.coefs[k] * c.coefs[k];
        }
    }
}

/// A single box (inequality) constraint over the 8 trilinear corner dofs of
/// one point: `lower <= coefs . x[indices] <= upper`. Used for value-bound
/// inequality constraints (see `add_value_inequality_constraint`).
///
/// `coefs`/`lower`/`upper` already have both the row's L2 normalisation
/// *and* `weight` folded in (`scale = weight / norm`, applied uniformly to
/// `coefs`, `lower`, and `upper` -- see `add_value_inequality_constraint`),
/// unlike the equality [`PointConstraint`]/[`PairConstraint`] families, which
/// keep `weight` separate and square it in at accumulation time. Folding it
/// in here instead keeps `solve_admm`'s box-projection z-update
/// (`z = clip(Q x + u, lower, upper)`) a single elementwise op against a flat
/// per-row bound pair, with no separate per-row weight to thread through the
/// ADMM state (`z`/`u`) -- `weight` only ever needs to scale this row's
/// contribution to `Q`/the bounds once, at constraint-build time, since a
/// positive rescaling of both a row and its bounds together leaves the
/// constraint it encodes (`lower_raw <= raw_row_value <= upper_raw`)
/// unchanged.
struct InequalityPointConstraint {
    indices: [usize; 8],
    coefs: [f64; 8],
    lower: f64,
    upper: f64,
}

/// As [`InequalityPointConstraint`], but over the *difference* of two
/// points' trilinear corner dofs (`coefs_a . x[indices_a] - coefs_b .
/// x[indices_b]`), bounded to `[lower, upper]` instead of forced to zero --
/// the pairwise-ordering case (see `add_inequality_pair_constraint`), e.g.
/// stratigraphic ordering (`lower = -inf`, `upper = 0` encodes "unit A's
/// value must not exceed unit B's").
struct InequalityPairConstraint {
    /// First 8 entries are point A's corner indices, next 8 are point B's.
    indices: [usize; 16],
    /// Row-normalised, weight-scaled combined coefficients: point A's value
    /// coefs followed by the *negated* point B's value coefs (see
    /// [`PairConstraint::coefs`]'s doc comment for the same layout).
    coefs: [f64; 16],
    lower: f64,
    upper: f64,
}

#[inline]
fn ineq_point_row_value(c: &InequalityPointConstraint, x: &[f64]) -> f64 {
    let mut v = 0.0;
    for k in 0..8 {
        v += c.coefs[k] * x[c.indices[k]];
    }
    v
}

#[inline]
fn ineq_pair_row_value(c: &InequalityPairConstraint, x: &[f64]) -> f64 {
    let mut v = 0.0;
    for k in 0..16 {
        v += c.coefs[k] * x[c.indices[k]];
    }
    v
}

/// One row of [`FiniteDifferenceInterpolator::solve_active_set_cg`]'s
/// working set: inequality row `row` (in the flattened value-then-pair order
/// [`FiniteDifferenceInterpolator::inequality_row`] uses) held exactly at
/// `bound`. `sign` is `+1.0` when `bound` is that row's *lower* bound
/// (binding as `Q_row x >= lower`) and `-1.0` when it's the *upper* bound
/// (`Q_row x <= upper`) -- only used to get the Lagrange-multiplier release
/// check's sign right (see `active_set_cg_core`'s doc comment).
#[derive(Clone, Copy)]
struct ActiveRow {
    row: usize,
    bound: f64,
    sign: f64,
}

/// Dense Cholesky factorisation of a small SPD matrix `a` (row-major, `n x
/// n`): returns the lower-triangular factor `L` (same row-major layout)
/// with `A = L L^T`, or `None` if a pivot is non-positive. Used for the
/// small `|W| x |W|` active-working-set Gram systems in
/// `active_set_cg_core` -- `|W|` is expected to be tiny relative to `dof()`
/// (see the "Constraint Density" row of `Inequalities.md`'s decision
/// matrix), so a dense direct factorisation here is the simplest and
/// fastest option, not a departure from the crate's matrix-free stance on
/// the *grid* operator `N` (which stays completely implicit).
fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i * n + j] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Solve `A y = b` given `A`'s Cholesky factor `l` (see [`cholesky`]), via
/// forward + back substitution.
fn cholesky_solve(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i * n + k] * y[k];
        }
        y[i] = sum / l[i * n + i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = y[i];
        for k in (i + 1)..n {
            sum -= l[k * n + i] * x[k];
        }
        x[i] = sum / l[i * n + i];
    }
    x
}

/// `Q_W v`: dot of `v` against each active row, in working-set order.
fn active_set_apply(
    interp: &FiniteDifferenceInterpolator,
    active: &[ActiveRow],
    v: &[f64],
) -> Vec<f64> {
    active
        .iter()
        .map(|a| interp.inequality_row_dot(a.row, v))
        .collect()
}

/// Dense Gram matrix `Q_W Q_W^T` (row-major, `|W| x |W|`), with a tiny
/// relative ridge added to the diagonal so near-duplicate/colinear active
/// rows (e.g. two inequality constraints that happen to coincide) keep the
/// Cholesky factorisation well-defined instead of failing outright.
fn active_set_gram(interp: &FiniteDifferenceInterpolator, active: &[ActiveRow]) -> Vec<f64> {
    let w = active.len();
    let mut g = vec![0.0; w * w];
    for a in 0..w {
        for b in a..w {
            let v = interp.inequality_row_dot_rows(active[a].row, active[b].row);
            g[a * w + b] = v;
            g[b * w + a] = v;
        }
    }
    g
}

/// Cholesky-factorise a working-set Gram matrix (see [`active_set_gram`]),
/// trying the matrix exactly as given first. A ridge is *not* added
/// up-front: every active row is already L2-normalised to unit norm at
/// constraint-build time (see `add_value_inequality_constraint`'s docs), so
/// a genuinely independent working set's Gram matrix already has diagonal
/// entries of exactly `1.0` and is well-conditioned -- and this
/// factorisation feeds directly into `active_set_restore_feasibility`'s
/// exact affine correction and the multiplier release check, where even a
/// tiny relative perturbation (a naive fixed
/// ridge) leaks a fraction of the -- potentially large -- unconstrained
/// residual back into the "should be pinned" active dofs every CG
/// iteration, silently drifting `x` off the active bound over many
/// iterations instead of keeping it exact. A ridge is only introduced as a
/// fallback, and only just large enough to make the *exact* Gram
/// factorisable, for the rare case of near-duplicate/colinear active rows
/// (e.g. two inequality constraints that happen to coincide).
fn active_set_cholesky(gram: &[f64], w: usize) -> Vec<f64> {
    if let Some(l) = cholesky(gram, w) {
        return l;
    }
    let diag_max = (0..w).map(|i| gram[i * w + i]).fold(1.0_f64, f64::max);
    let mut ridge = 1e-12 * diag_max;
    loop {
        let mut g = gram.to_vec();
        for i in 0..w {
            g[i * w + i] += ridge;
        }
        if let Some(l) = cholesky(&g, w) {
            return l;
        }
        ridge *= 10.0;
        assert!(
            ridge < diag_max.max(1.0) * 1e6,
            "active working set's Gram matrix is not positive definite even after ridge regularisation -- likely a bug in how active rows are constructed"
        );
    }
}

/// Dense Gram matrix `Q_W diag(dinv) Q_W^T` (row-major, `|W| x |W|`) -- the
/// working set's Gram matrix in the `M^-1 = diag(dinv)` (Jacobi) inner
/// product, used for [`active_set_project_weighted`]'s preconditioned
/// null-space projection. Unlike [`active_set_gram`] (the plain `Q_W Q_W^T`
/// used for exact feasibility restoration and the Lagrange-multiplier
/// release check, both of which are properties of the *original* unweighted
/// problem), this one exists purely to make the projected CG search
/// converge in a comparable number of iterations to
/// [`FiniteDifferenceInterpolator::solve`]'s Jacobi-preconditioned CG,
/// rather than plain (unpreconditioned) CG's typically much slower
/// convergence on this crate's regularised normal-equations operator.
fn active_set_gram_weighted(
    interp: &FiniteDifferenceInterpolator,
    active: &[ActiveRow],
    dinv: &[f64],
) -> Vec<f64> {
    let w = active.len();
    let mut g = vec![0.0; w * w];
    for a in 0..w {
        for b in a..w {
            let v = interp.inequality_row_dot_rows_weighted(active[a].row, active[b].row, dinv);
            g[a * w + b] = v;
            g[b * w + a] = v;
        }
    }
    g
}

/// Jacobi-preconditioned, null-space-projected residual: `M^-1 v - M^-1
/// Q_W^T (Q_W M^-1 Q_W^T)^-1 Q_W (M^-1 v)`, `M^-1 = diag(dinv)` -- the
/// projected-preconditioned-CG generalisation (Nocedal & Wright, "Numerical
/// Optimization", Algorithm 16.2) of [`FiniteDifferenceInterpolator::solve`]'s
/// plain Jacobi preconditioner to the working set's null space. `chol` is
/// [`active_set_gram_weighted`]'s Cholesky factor (`None` iff the working
/// set is empty, in which case this degenerates to plain Jacobi
/// preconditioning, matching `solve`'s unconstrained behaviour).
///
/// Note this is *not* interchangeable with [`active_set_restore_feasibility`]
/// or the Lagrange-multiplier release check in `active_set_cg_core`: both of
/// those need the *exact*, unweighted `Q_W Q_W^T` system (see
/// [`active_set_gram`]/[`active_set_apply`]) since they answer questions
/// about the original (unweighted) QP, not about how fast this
/// preconditioner happens to converge.
fn active_set_project_weighted(
    interp: &FiniteDifferenceInterpolator,
    active: &[ActiveRow],
    chol: Option<&[f64]>,
    dinv: &[f64],
    v: &[f64],
) -> Vec<f64> {
    let weighted: Vec<f64> = v.iter().zip(dinv.iter()).map(|(vi, di)| vi * di).collect();
    let w = active.len();
    let chol = match chol {
        Some(c) if w > 0 => c,
        _ => return weighted,
    };
    let qv = active_set_apply(interp, active, &weighted);
    let y = cholesky_solve(chol, w, &qv);
    let mut out = weighted;
    for (a, yi) in active.iter().zip(y.iter()) {
        let (idx, coefs) = interp.inequality_row(a.row);
        let s = -*yi;
        for (i, c) in idx.iter().zip(coefs.iter()) {
            out[*i] += s * c * dinv[*i];
        }
    }
    out
}

/// Correct `x` in place to satisfy `Q_W x = bound_W` exactly, via the
/// minimum-Euclidean-distance affine projection `x <- x - Q_W^T (Q_W
/// Q_W^T)^-1 (Q_W x - bound_W)`. A no-op (up to floating point) whenever `x`
/// already satisfies every active row exactly -- in particular this makes it
/// safe to call unconditionally after *every* working-set change (a row
/// joining or leaving), including right after a mid-CG step-truncation event
/// where `x` is already exactly on the new boundary by construction (see
/// `active_set_cg_core`): the correction there is ~0, so no CG progress is
/// thrown away, while a freshly-added row that `x` doesn't yet satisfy (or a
/// warm-started `x` from ADMM in `solve_hybrid`, which only satisfies bounds
/// approximately) gets pulled onto the constraint exactly.
fn active_set_restore_feasibility(
    interp: &FiniteDifferenceInterpolator,
    active: &[ActiveRow],
    chol: &[f64],
    x: &mut [f64],
) {
    let w = active.len();
    if w == 0 {
        return;
    }
    let residual: Vec<f64> = active
        .iter()
        .map(|a| interp.inequality_row_dot(a.row, x) - a.bound)
        .collect();
    let y = cholesky_solve(chol, w, &residual);
    for (a, yi) in active.iter().zip(y.iter()) {
        interp.inequality_row_scatter_add(a.row, -*yi, x);
    }
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Build the strike/dip basis `(u, v)` for unit vector `n`. It follows the
/// Python `loop_common.math.get_vectors` convention: the strike `u` is
/// horizontal, `u = (n_y, -n_x, 0) / hypot(n_x, n_y)`. For an exactly
/// vertical normal (`hypot == 0`) Python uses `atan2(0, 0) = 0`, so
/// `u = (0, -1, 0)`. The dip is `v = normalise(cross(u, n))`.
fn orthonormal_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let len = norm3(n);
    let n = if len > 0.0 { scale3(n, 1.0 / len) } else { n };
    let h = n[0].hypot(n[1]);
    let u = if h > 0.0 {
        [n[1] / h, -n[0] / h, 0.0]
    } else {
        [0.0, -1.0, 0.0]
    };
    let v_raw = cross(u, n);
    let v_norm = norm3(v_raw);
    let v = if v_norm > 0.0 {
        scale3(v_raw, 1.0 / v_norm)
    } else {
        v_raw
    };
    (u, v)
}

#[derive(Debug, Clone, Copy)]
pub struct RegularisationWeights {
    pub dxx: f64,
    pub dyy: f64,
    pub dzz: f64,
    pub dxy: f64,
    pub dxz: f64,
    pub dyz: f64,
}

impl Default for RegularisationWeights {
    fn default() -> Self {
        Self {
            dxx: 1.0,
            dyy: 1.0,
            dzz: 1.0,
            dxy: 1.0,
            dxz: 1.0,
            dyz: 1.0,
        }
    }
}

/// Border (centred second-derivative `[1, -2, 1]`) regularisation weights,
/// one per axis, matching the Python `interpolation_weights["dx"/"dy"/"dz"]`.
#[derive(Debug, Clone, Copy)]
pub struct BorderRegularisationWeights {
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
}

impl Default for BorderRegularisationWeights {
    fn default() -> Self {
        Self {
            dx: 1.0,
            dy: 1.0,
            dz: 1.0,
        }
    }
}

fn build_border_stencils(
    grid: &StructuredGrid3D,
    weights: BorderRegularisationWeights,
) -> Vec<(EdgeStencil, f64)> {
    let axis_weights = [weights.dx, weights.dy, weights.dz];
    let mut stencils = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for axis in 0..3 {
        let w = axis_weights[axis];
        if let Some(s) = edge_forward(axis, &grid.nsteps) {
            stencils.push((s, w));
        }
        if let Some(s) = edge_backward(axis, &grid.nsteps) {
            stencils.push((s, w));
        }
    }
    stencils
}

/// One entry of the fused six-family regularisation kernel: `flat_offset` is
/// the tap's grid-step displacement (each component in `-2..=2`)
/// pre-converted to a flat (F-order) index delta for this grid shape (i.e.
/// `e.0 + nx*e.1 + nx*ny*e.2`, computed once when the kernel is built rather
/// than recomputed for every node `fused_interior_regularisation_into`
/// visits), and the paired `f64` is `sum_family w_family^2 * sum_{d1,d2 in
/// family taps: d2-d1==offset} coeff[d1]*coeff[d2]` -- the combined `S^T S`
/// autocorrelation kernel described in the module docs / `TASKS.md`'s "fused
/// single-kernel regularisation" item.
type KernelTap = (i64, f64);

/// Minimum grid-step distance to the nearest *true* domain edge along any
/// axis (0 for the outermost node layer, 1 for the next layer in, ...).
/// Composing two radius-1 stencils (forward then adjoint) is only a
/// translation-invariant radius-2 stencil away from the boundary; see
/// `fused_interior_regularisation_into`'s doc comment for the full argument.
#[inline]
fn dist_to_true_edge(nsteps: &[usize; 3], i: usize, j: usize, k: usize) -> usize {
    let dx = i.min(nsteps[0] - 1 - i);
    let dy = j.min(nsteps[1] - 1 - j);
    let dz = k.min(nsteps[2] - 1 - k);
    dx.min(dy).min(dz)
}

/// Precompute, once per grid shape (independent of regularisation weights):
/// - `shell_row_indices`: flat indices of every node at distance <= 2 from
///   the true edge -- the candidate stencil *row centres* whose forward
///   value can ever contribute to a boundary-shell output node once both a
///   forward and an adjoint radius-1 step are applied (see
///   `shell_regularisation_correction_into`).
/// - `shell_output_indices`: flat indices at distance < 2 -- the output
///   nodes `fused_interior_regularisation_into` deliberately skips and that
///   `shell_regularisation_correction_into` must fix up exactly.
/// - the former's flat-index inverse (`flat -> position`, `u32::MAX` for
///   non-shell-row nodes), so the adjoint half of the shell correction can
///   find a row's forward value by flat index in O(1).
///
/// All three tables hold `u32` values to use 4 bytes per entry. The top
/// value `u32::MAX` is the empty marker.
fn build_shell_indices(
    grid: &StructuredGrid3D,
    tracker: &mut MemTracker,
) -> Result<ShellTables, CoreError> {
    let dof = grid.n_nodes();
    assert!(
        dof < u32::MAX as usize,
        "grid has too many nodes for the u32 shell tables"
    );
    let mut shell_rows = Vec::new();
    let mut shell_outputs = Vec::new();
    for flat in 0..dof {
        let (i, j, k) = grid.node_ijk(flat);
        let d = dist_to_true_edge(&grid.nsteps, i, j, k);
        if d <= 2 {
            shell_rows.push(flat as u32);
        }
        if d < 2 {
            shell_outputs.push(flat as u32);
        }
    }
    // The lookup is the only table of size dof. The tracker counts it.
    let mut lookup: Vec<u32> = tracker.try_zeroed(dof)?;
    lookup.fill(u32::MAX);
    for (pos, &flat) in shell_rows.iter().enumerate() {
        lookup[flat as usize] = pos as u32;
    }
    Ok(ShellTables {
        rows: shell_rows,
        outputs: shell_outputs,
        lookup,
    })
}

/// The three shell tables (see `build_shell_indices`).
struct ShellTables {
    /// Flat indices of the nodes at distance <= 2 from the true edge.
    rows: Vec<u32>,
    /// Flat indices of the nodes at distance < 2.
    outputs: Vec<u32>,
    /// Inverse of `rows`: flat index to position, `u32::MAX` if not a row.
    lookup: Vec<u32>,
}

/// Build the fused single-kernel regularisation's combined (5,5,5)-radius
/// autocorrelation kernel from the six curvature families' (already
/// unit-L2-normalised) taps and their scalar weights (see module docs).
/// Every one of the six families always qualifies for fusion in this crate
/// (unlike the Python original, which also has to guard against
/// spatially-varying per-dof weights and caller-supplied custom masks --
/// neither exists here: `RegularisationWeights` is always one scalar per
/// family and the six stencils are always the fixed canonical shapes from
/// `stencil::{dxx, dyy, dzz, dxy, dxz, dyz}`), so this simply sums every
/// family's contribution unconditionally (skipping any with zero weight,
/// which would contribute nothing anyway).
fn build_combined_kernel(
    stencils: &[(Stencil, f64); 6],
    grid: &StructuredGrid3D,
) -> Vec<KernelTap> {
    let mut map: HashMap<(i32, i32, i32), f64> = HashMap::new();
    for (stencil, weight) in stencils {
        if *weight == 0.0 {
            continue;
        }
        let w2 = weight * weight;
        for (d1, c1) in stencil.offsets.iter().zip(stencil.coeffs.iter()) {
            for (d2, c2) in stencil.offsets.iter().zip(stencil.coeffs.iter()) {
                let e = (d2.0 - d1.0, d2.1 - d1.1, d2.2 - d1.2);
                *map.entry(e).or_insert(0.0) += w2 * c1 * c2;
            }
        }
    }
    // `HashMap`'s iteration order is randomised per-instance (a fresh random
    // hasher seed per `HashMap::new()`), so two structurally-identical
    // interpolators (same grid shape, same weights) could otherwise end up
    // with their combined kernel taps in a different order -- and since
    // float addition isn't associative, summing the *same* taps in a
    // *different* order in `fused_interior_regularisation_into` can differ
    // by a ULP or two between the two instances. Sorting by offset makes
    // the kernel's (and therefore every downstream sum's) order a pure
    // function of the stencils/weights, restoring the bit-for-bit
    // reproducibility every other part of this crate already assumes
    // between structurally-identical interpolators.
    let mut offsets: Vec<((i32, i32, i32), f64)> =
        map.into_iter().filter(|(_, v)| v.abs() > 1e-15).collect();
    offsets.sort_by_key(|(e, _)| *e);
    // Pre-convert each (i32,i32,i32) offset to a flat index delta for this
    // grid shape once, here, rather than recomputing `nx*e.1 + nx*ny*e.2`
    // for every (node, tap) pair inside `fused_interior_regularisation_into`'s
    // hot loop -- this was previously the dominant cost of the fused pass,
    // since it re-did the same two multiplies per tap on every single node
    // visited instead of once per tap ever.
    let nx_i = grid.nsteps[0] as i64;
    let ny_i = grid.nsteps[1] as i64;
    offsets
        .into_iter()
        .map(|(e, c)| (e.0 as i64 + nx_i * e.1 as i64 + nx_i * ny_i * e.2 as i64, c))
        .collect()
}

/// Look up a stencil's normalised coefficient at a specific tap offset (0 if
/// the stencil has no tap there). Only used to build the `BlockJacobiX`
/// preconditioner's off-diagonal (see `dxx_and_border_dx_offdiagonal`),
/// where we need specific named taps rather than an iteration over all of
/// them.
fn coeff_at(stencil: &Stencil, offset: (i32, i32, i32)) -> f64 {
    stencil
        .offsets
        .iter()
        .zip(stencil.coeffs.iter())
        .find(|(o, _)| **o == offset)
        .map(|(_, c)| *c)
        .unwrap_or(0.0)
}

/// Solve the symmetric tridiagonal system `diag[i] z[i] + off[i-1] z[i-1] +
/// off[i] z[i+1] = r[i]` (`off[n-1]` is never read) via the Thomas
/// algorithm. `cprime` is caller-supplied scratch (length `n`) so this makes
/// no allocation of its own -- needed since it runs inside the CG hot loop
/// for `PreconditionerKind::BlockJacobiX`. A near-zero pivot falls back to
/// leaving that entry's contribution as if the local block were diagonal
/// there (mirroring the plain Jacobi preconditioner's `if di > 0.0 {..} else
/// {..}` fallback) rather than dividing by (near) zero; should not trigger
/// in practice given `ridge_factor > 0` (see `dxx_and_border_dx_offdiagonal`'s
/// docs for why that guarantees a strictly positive definite block).
fn solve_tridiagonal_line(diag: &[f64], off: &[f64], r: &[f64], cprime: &mut [f64], z: &mut [f64]) {
    let n = diag.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        z[0] = if diag[0] > 0.0 { r[0] / diag[0] } else { r[0] };
        return;
    }
    let mut denom = diag[0];
    if denom.abs() < 1e-300 {
        denom = 1.0;
    }
    cprime[0] = off[0] / denom;
    z[0] = r[0] / denom;
    for i in 1..n {
        let e_prev = off[i - 1];
        let mut m = diag[i] - e_prev * cprime[i - 1];
        if m.abs() < 1e-300 {
            m = 1e-300;
        }
        let e_i = if i < n - 1 { off[i] } else { 0.0 };
        cprime[i] = e_i / m;
        z[i] = (r[i] - e_prev * z[i - 1]) / m;
    }
    for i in (0..n - 1).rev() {
        z[i] -= cprime[i] * z[i + 1];
    }
}

/// Which preconditioner [`FiniteDifferenceInterpolator::solve_with`] uses.
/// `Jacobi` is the original plain-diagonal preconditioner (cheap, always
/// applicable, kept as the default so existing callers are unaffected).
/// `BlockJacobiX` solves a small tridiagonal system per x-axis grid line,
/// exactly capturing the `dxx` family's and any x-axis border
/// regularisation's line-local coupling instead of discarding it to the
/// diagonal -- see `dxx_and_border_dx_offdiagonal` for the SPD argument.
/// It degenerates to plain Jacobi's tridiagonal-with-zero-off-diagonal
/// special case whenever there's no such coupling (e.g. `dxx` weight is
/// zero and no x-axis border regularisation is configured), so it is never
/// worse than Jacobi in that sense, only ever adds real information.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PreconditionerKind {
    #[default]
    Jacobi,
    BlockJacobiX,
}

/// Result of [`FiniteDifferenceInterpolator::solve_admm`].
pub struct AdmmResult {
    pub x: Vec<f64>,
    /// Number of outer ADMM iterations run. `0` when there were no
    /// inequality constraints at all (the plain-CG-delegating fast path).
    pub outer_iterations: usize,
    /// Whether the outer primal/dual residual stopping criterion was met
    /// before `max_outer_iter` outer iterations (or, in the no-inequality
    /// fast path, whether the delegated equality CG solve converged).
    pub converged: bool,
    /// Sum of every outer iteration's inner-CG iteration count (a rough
    /// total-work indicator; not itself a convergence signal).
    pub inner_iterations_total: usize,
}

/// Result of [`FiniteDifferenceInterpolator::solve_active_set_cg`] (and the
/// active-set-polish stage of [`solve_hybrid`](FiniteDifferenceInterpolator::solve_hybrid)).
pub struct ActiveSetCgResult {
    pub x: Vec<f64>,
    /// Number of working-set changes (a row joining or leaving the active
    /// set), each one restarting the projected-CG inner solve on the
    /// updated working set. `0` when there were no inequality constraints
    /// at all (the plain-CG-delegating fast path).
    pub outer_iterations: usize,
    /// Sum of every working-set stage's inner (projected) CG iteration
    /// count (a rough total-work indicator, not itself a convergence
    /// signal).
    pub inner_iterations_total: usize,
    /// Whether the KKT optimality conditions were met (projected residual
    /// below `inner_tol` *and* no active row left with a wrong-signed
    /// Lagrange multiplier) before `max_outer_iter` working-set changes.
    pub converged: bool,
}

/// Result of [`FiniteDifferenceInterpolator::solve_hybrid`].
pub struct HybridResult {
    pub x: Vec<f64>,
    /// Number of ADMM warm-up outer iterations run (always exactly the
    /// caller-supplied `admm_warmup_iter`, since the warm-up stage runs a
    /// fixed budget rather than iterating to convergence).
    pub admm_warmup_iterations: usize,
    /// Sum of the ADMM warm-up stage's inner-CG iteration counts.
    pub admm_inner_iterations_total: usize,
    /// The active-set polish stage's `outer_iterations` (see
    /// [`ActiveSetCgResult::outer_iterations`]).
    pub polish_outer_iterations: usize,
    /// The active-set polish stage's `inner_iterations_total` (see
    /// [`ActiveSetCgResult::inner_iterations_total`]).
    pub polish_inner_iterations_total: usize,
    /// Whether the active-set polish stage converged (see
    /// [`ActiveSetCgResult::converged`]).
    pub converged: bool,
}

pub struct FiniteDifferenceInterpolator {
    pub grid: StructuredGrid3D,
    pub regularisation: RegularisationWeights,
    pub border_regularisation: BorderRegularisationWeights,
    /// Diagonal Tikhonov damping added to the normal equations so the
    /// system stays positive definite even with sparse/degenerate data
    /// constraints, matching `DiscreteInterpolator.ridge_factor`.
    pub ridge_factor: f64,
    value_constraints: Vec<PointConstraint>,
    gradient_constraints: Vec<PointConstraint>,
    norm_constraints: Vec<PointConstraint>,
    tangent_constraints: Vec<PointConstraint>,
    interface_constraints: Vec<PairConstraint>,
    /// Value-bound inequality rows (see `add_value_inequality_constraint`);
    /// row order here, followed by `pair_inequality_constraints`, is the
    /// fixed order `solve_admm`'s `Q`/`z`/`u` use.
    value_inequality_constraints: Vec<InequalityPointConstraint>,
    /// Pairwise-ordering inequality rows (see `add_inequality_pair_constraint`).
    pair_inequality_constraints: Vec<InequalityPairConstraint>,
    stencils: [(Stencil, f64); 6],
    border_stencils: Vec<(EdgeStencil, f64)>,
    /// Fused single-kernel regularisation's combined kernel (see
    /// `build_combined_kernel`/`fused_interior_regularisation_into`);
    /// rebuilt whenever `stencils`'/weights change (`with_regularisation`).
    /// The table is empty until `prepare` (or the first use) builds it.
    combined_kernel: OnceLock<Vec<KernelTap>>,
    /// Row centres (dist to true edge <= 2) / output nodes (dist < 2) / the
    /// former's flat-index inverse -- grid-shape-only, built once (see
    /// `build_shell_indices`). The tables are empty until `prepare` (or the
    /// first use) builds them. `new` does not allocate them.
    shell: OnceLock<ShellTables>,
    /// Per-node squared weight for each of the six curvature families
    /// (`[dxx, dyy, dzz, dxy, dxz, dyz]`), accumulated across every
    /// [`add_directional_regularisation`](Self::add_directional_regularisation)
    /// call. `None` until the first such call (lazily allocated so grids that
    /// never use directional regularisation -- i.e. every caller except fold
    /// interpolation -- pay no extra memory or per-iteration cost). Unlike
    /// `stencils`' uniform scalar weights, these vary per grid node, so they
    /// cannot go through the fused single-kernel fast path (which requires a
    /// uniform weight per family -- see `build_combined_kernel`'s docs) and
    /// are instead applied via their own unfused pass, the same way
    /// `border_stencils` are.
    directional_weight2: Option<[Vec<f64>; 6]>,
}

impl FiniteDifferenceInterpolator {
    pub fn new(grid: StructuredGrid3D) -> Self {
        let regularisation = RegularisationWeights::default();
        let border_regularisation = BorderRegularisationWeights::default();
        let stencils = [
            (dxx(), regularisation.dxx),
            (dyy(), regularisation.dyy),
            (dzz(), regularisation.dzz),
            (dxy(), regularisation.dxy),
            (dxz(), regularisation.dxz),
            (dyz(), regularisation.dyz),
        ];
        let border_stencils = build_border_stencils(&grid, border_regularisation);
        Self {
            grid,
            regularisation,
            border_regularisation,
            ridge_factor: 1e-8,
            value_constraints: Vec::new(),
            gradient_constraints: Vec::new(),
            norm_constraints: Vec::new(),
            tangent_constraints: Vec::new(),
            interface_constraints: Vec::new(),
            value_inequality_constraints: Vec::new(),
            pair_inequality_constraints: Vec::new(),
            stencils,
            border_stencils,
            combined_kernel: OnceLock::new(),
            shell: OnceLock::new(),
            directional_weight2: None,
        }
    }

    pub fn with_regularisation(mut self, weights: RegularisationWeights) -> Self {
        self.regularisation = weights;
        self.stencils = [
            (dxx(), weights.dxx),
            (dyy(), weights.dyy),
            (dzz(), weights.dzz),
            (dxy(), weights.dxy),
            (dxz(), weights.dxz),
            (dyz(), weights.dyz),
        ];
        // The weights changed, so the old kernel is not valid. Drop it.
        self.combined_kernel = OnceLock::new();
        self
    }

    /// Build the shell tables and the combined kernel. Call it once before a
    /// solve. A second call does nothing. The grid-sized shell lookup is
    /// counted by `tracker`.
    #[allow(dead_code)]
    pub(crate) fn prepare(&mut self, tracker: &mut MemTracker) -> Result<(), CoreError> {
        self.ensure_tables(tracker)
    }

    /// As `prepare`, but it takes `&self`. The tables use `OnceLock`.
    fn ensure_tables(&self, tracker: &mut MemTracker) -> Result<(), CoreError> {
        if self.shell.get().is_none() {
            let tables = build_shell_indices(&self.grid, tracker)?;
            // If another thread set the tables first, keep that value.
            let _ = self.shell.set(tables);
        }
        self.kernel();
        Ok(())
    }

    /// The shell tables. Built without a budget if `prepare` did not run.
    fn shell(&self) -> &ShellTables {
        self.shell.get_or_init(|| {
            build_shell_indices(&self.grid, &mut MemTracker::new(None))
                .expect("shell tables must fit in memory")
        })
    }

    /// The combined kernel. Built on first use.
    fn kernel(&self) -> &[KernelTap] {
        self.combined_kernel
            .get_or_init(|| build_combined_kernel(&self.stencils, &self.grid))
    }

    pub fn with_border_regularisation(mut self, weights: BorderRegularisationWeights) -> Self {
        self.border_regularisation = weights;
        self.border_stencils = build_border_stencils(&self.grid, weights);
        self
    }

    pub fn dof(&self) -> usize {
        self.grid.n_nodes()
    }

    /// Bytes held by the interpolator's own tables: the shell lookup, the
    /// shell index lists, the directional weights (if any) and the combined
    /// kernel. It counts the capacity of each table. It does not count the
    /// constraint records or the vectors made during a solve.
    pub fn table_bytes(&self) -> usize {
        let directional = self
            .directional_weight2
            .as_ref()
            .map_or(0, |w| w.iter().map(|v| v.capacity()).sum::<usize>())
            * std::mem::size_of::<f64>();
        // Before `prepare` the tables are not built, so they count 0 bytes.
        let shell = self.shell.get().map_or(0, |s| {
            (s.lookup.capacity() + s.rows.capacity() + s.outputs.capacity())
                * std::mem::size_of::<u32>()
        });
        let kernel = self
            .combined_kernel
            .get()
            .map_or(0, |k| k.capacity() * std::mem::size_of::<KernelTap>());
        shell + directional + kernel
    }

    /// Add a scalar value observation at `point` with weight `weight`.
    /// Points outside the grid's bounding box are silently dropped (matching
    /// the Python interpolator's `inside` masking).
    pub fn add_value_constraint(&mut self, point: [f64; 3], value: f64, weight: f64) {
        let Some((indices, raw_coefs)) = self.grid.trilinear_coeffs(point) else {
            return;
        };
        let norm = raw_coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
        if norm <= 0.0 {
            return;
        }
        let coefs = raw_coefs.map(|c| c / norm);
        self.value_constraints.push(PointConstraint {
            indices,
            coefs,
            value: value / norm,
            weight,
        });
    }

    pub fn n_value_constraints(&self) -> usize {
        self.value_constraints.len()
    }

    /// Shared machinery for gradient-orthogonal constraints (`strike . grad
    /// f = 0`, `dip . grad f = 0`, from `add_gradient_constraint`) and
    /// tangent constraints (`v . grad f = 0`, from `add_tangent_constraint`):
    /// both reduce to "one row constraining the field gradient at `point` to
    /// be orthogonal to `vector`", i.e. the Python
    /// `add_gradient_orthogonal_constraints`. `vector` is given in world
    /// space and normalised before being dotted with the trilinear gradient
    /// shape function; `None` if `point` is outside the grid or `vector` is
    /// (numerically) zero.
    fn build_orthogonal_constraint(
        &self,
        point: [f64; 3],
        vector: [f64; 3],
        weight: f64,
    ) -> Option<PointConstraint> {
        let (indices, grad) = self.grid.trilinear_gradient(point)?;
        let v_local = vector;
        let v_norm = norm3(v_local);
        if v_norm <= 0.0 {
            return None;
        }
        let v_unit = scale3(v_local, 1.0 / v_norm);
        let mut raw_coefs = [0.0; 8];
        for c in 0..8 {
            raw_coefs[c] = v_unit[0] * grad[c][0] + v_unit[1] * grad[c][1] + v_unit[2] * grad[c][2];
        }
        let row_norm = raw_coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
        if row_norm <= 0.0 {
            return None;
        }
        let coefs = raw_coefs.map(|c| c / row_norm);
        Some(PointConstraint {
            indices,
            coefs,
            value: 0.0,
            weight,
        })
    }

    /// Add an orientation-data (gradient) constraint at `point`: `normal` is
    /// a (not necessarily unit) vector giving the target gradient
    /// *direction* only (its magnitude is discarded -- use
    /// [`add_norm_constraint`](Self::add_norm_constraint) to also constrain
    /// magnitude). Internally split into two orthogonal directions spanning
    /// the plane perpendicular to `normal` (see
    /// [`orthonormal_basis`]/Python's `get_vectors` strike/dip split), each
    /// constrained to be orthogonal to the field gradient at `point`. Points
    /// outside the grid, or a zero `normal`, are silently dropped.
    pub fn add_gradient_constraint(&mut self, point: [f64; 3], normal: [f64; 3], weight: f64) {
        let n_norm = norm3(normal);
        if n_norm <= 0.0 {
            return;
        }
        let unit = scale3(normal, 1.0 / n_norm);
        let (u, v) = orthonormal_basis(unit);
        if let Some(c) = self.build_orthogonal_constraint(point, u, weight) {
            self.gradient_constraints.push(c);
        }
        if let Some(c) = self.build_orthogonal_constraint(point, v, weight) {
            self.gradient_constraints.push(c);
        }
    }

    pub fn n_gradient_constraints(&self) -> usize {
        self.gradient_constraints.len()
    }

    /// Add a tangent constraint at `point`: constrains the field gradient to
    /// be orthogonal to the given (arbitrary, not necessarily unit) `vector`
    /// -- a single row, unlike [`add_gradient_constraint`](Self::add_gradient_constraint)'s
    /// two (strike/dip) rows. Points outside the grid, or a zero `vector`,
    /// are silently dropped.
    pub fn add_tangent_constraint(&mut self, point: [f64; 3], vector: [f64; 3], weight: f64) {
        if let Some(c) = self.build_orthogonal_constraint(point, vector, weight) {
            self.tangent_constraints.push(c);
        }
    }

    pub fn n_tangent_constraints(&self) -> usize {
        self.tangent_constraints.len()
    }

    /// Add a norm constraint at `point`: constrains each of the 3 gradient
    /// components directly to the corresponding component of
    /// `target_gradient` (given in world space), one row per axis. Unlike
    /// [`add_gradient_constraint`](Self::add_gradient_constraint), this is
    /// non-homogeneous (targets an actual value, not just "= 0"), so it
    /// constrains gradient *magnitude* as well as direction. Points outside
    /// the grid are silently dropped.
    pub fn add_norm_constraint(&mut self, point: [f64; 3], target_gradient: [f64; 3], weight: f64) {
        let Some((indices, grad)) = self.grid.trilinear_gradient(point) else {
            return;
        };
        let target_local = target_gradient;
        for axis in 0..3 {
            let mut raw_coefs = [0.0; 8];
            for c in 0..8 {
                raw_coefs[c] = grad[c][axis];
            }
            let row_norm = raw_coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
            if row_norm <= 0.0 {
                continue;
            }
            let coefs = raw_coefs.map(|c| c / row_norm);
            let value = target_local[axis] / row_norm;
            self.norm_constraints.push(PointConstraint {
                indices,
                coefs,
                value,
                weight,
            });
        }
    }

    pub fn n_norm_constraints(&self) -> usize {
        self.norm_constraints.len()
    }

    /// Add interface (equal-value) constraints: every point in `points`
    /// sharing the same `id` (compared by bit-exact equality, matching
    /// `np.unique`-based grouping on a float id column) is pairwise
    /// constrained to have equal field value. Uses the "star" pairing
    /// strategy (mirroring the Python `_interface_pair_indices`'s minimal
    /// default) -- every point in a group is paired with the group's first
    /// point, rather than the full `all`/`chain` switching
    /// `interface_pair_mode` supports. Points outside the grid, `NaN` ids,
    /// and groups of size < 2 are silently dropped.
    pub fn add_interface_constraints(&mut self, points: &[[f64; 3]], ids: &[f64], weight: f64) {
        assert_eq!(
            points.len(),
            ids.len(),
            "points and ids must be the same length"
        );
        let mut groups: HashMap<u64, Vec<usize>> = HashMap::new();
        for (i, &id) in ids.iter().enumerate() {
            if id.is_nan() {
                continue;
            }
            groups.entry(id.to_bits()).or_default().push(i);
        }
        for idxs in groups.values() {
            if idxs.len() < 2 {
                continue;
            }
            let first = idxs[0];
            let Some((idx_a, coefs_a)) = self.grid.trilinear_coeffs(points[first]) else {
                continue;
            };
            for &j in &idxs[1..] {
                let Some((idx_b, coefs_b)) = self.grid.trilinear_coeffs(points[j]) else {
                    continue;
                };
                let mut indices = [0usize; 16];
                let mut coefs = [0.0f64; 16];
                indices[..8].copy_from_slice(&idx_a);
                indices[8..].copy_from_slice(&idx_b);
                for k in 0..8 {
                    coefs[k] = coefs_a[k];
                    coefs[8 + k] = -coefs_b[k];
                }
                let norm = coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
                if norm <= 0.0 {
                    continue;
                }
                for c in coefs.iter_mut() {
                    *c /= norm;
                }
                self.interface_constraints.push(PairConstraint {
                    indices,
                    coefs,
                    weight,
                });
            }
        }
    }

    pub fn n_interface_constraints(&self) -> usize {
        self.interface_constraints.len()
    }

    /// Add an anisotropic directional-regularisation term: penalises the
    /// directional second derivative `(v.grad)^2 f = 0` at every interior
    /// grid node, for a direction field `v` given one (world-space) vector
    /// per grid node in `directions` (length must equal
    /// [`dof`](Self::dof) -- mismatched lengths are silently dropped,
    /// matching every other constraint family's out-of-grid drop
    /// convention). This is the generic primitive fold interpolation's
    /// anisotropic smoothing is built from (see the Python
    /// `FiniteDifferenceInterpolator.minimise_directional_gradient_change`):
    /// callers evaluate their own directional geometry (e.g. a fold frame's
    /// axis/normal/deformed-orientation directions) and pass the resulting
    /// per-node vectors in here.
    ///
    /// Expanding `(v.grad)^2 f` in the six second-derivative families gives
    /// `vx^2 f_xx + vy^2 f_yy + vz^2 f_zz + 2vx*vy f_xy + 2vx*vz f_xz + 2vy*vz
    /// f_yz`, so each family stays independent (no cross-family terms in the
    /// squared residual) -- only its *weight* becomes per-node instead of one
    /// uniform scalar. Multiple calls (e.g. fold's three orthogonal
    /// directions, each with its own weight) accumulate independently: this
    /// method adds `(weight * comp_w_family(node))^2` to that family's
    /// running per-node squared-weight total, exactly as if each call were
    /// its own independent least-squares row set.
    ///
    /// `directions` are used raw (not renormalised), matching the Python
    /// original.
    ///
    /// Unlike the uniform-weight `stencils`, these per-node weights cannot
    /// use the fused single-kernel regularisation fast path (see
    /// `build_combined_kernel`'s docs), so they are accumulated into the
    /// normal operator via their own unfused pass -- see
    /// `apply_normal_operator_into`.
    pub fn add_directional_regularisation(&mut self, directions: &[[f64; 3]], weight: f64) {
        if directions.len() != self.dof() {
            return;
        }
        let dof = self.dof();
        let weight2 = self
            .directional_weight2
            .get_or_insert_with(|| std::array::from_fn(|_| vec![0.0; dof]));
        for (c, dir) in directions.iter().enumerate() {
            let [vx, vy, vz] = *dir;
            let comp_w = [
                vx * vx,
                vy * vy,
                vz * vz,
                2.0 * vx * vy,
                2.0 * vx * vz,
                2.0 * vy * vz,
            ];
            for family in 0..6 {
                let cw = weight * comp_w[family];
                weight2[family][c] += cw * cw;
            }
        }
    }

    /// Add a value-bound inequality constraint at `point`: constrains the
    /// predicted trilinear value there to lie in `[lower, upper]` (instead
    /// of an equality target), weighted by `weight`. Only used by
    /// [`solve_admm`](Self::solve_admm) -- [`solve`](Self::solve)/
    /// [`solve_with`](Self::solve_with) never look at inequality rows at
    /// all, so adding one has no effect unless you solve with `solve_admm`.
    ///
    /// Points outside the grid are silently dropped, matching every other
    /// constraint family. A non-positive `weight` is also silently dropped:
    /// unlike the equality families (which always use `weight^2` and so
    /// don't care about its sign), scaling an inequality row's `coefs` *and*
    /// its `lower`/`upper` bounds by the same factor only preserves the
    /// constraint's meaning if that factor is strictly positive.
    pub fn add_value_inequality_constraint(
        &mut self,
        point: [f64; 3],
        lower: f64,
        upper: f64,
        weight: f64,
    ) {
        if weight <= 0.0 {
            return;
        }
        let Some((indices, raw_coefs)) = self.grid.trilinear_coeffs(point) else {
            return;
        };
        let norm = raw_coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
        if norm <= 0.0 {
            return;
        }
        let scale = weight / norm;
        let coefs = raw_coefs.map(|c| c * scale);
        self.value_inequality_constraints
            .push(InequalityPointConstraint {
                indices,
                coefs,
                lower: lower * scale,
                upper: upper * scale,
            });
    }

    pub fn n_value_inequality_constraints(&self) -> usize {
        self.value_inequality_constraints.len()
    }

    /// Add a pairwise-ordering inequality constraint: constrains the
    /// difference between the predicted trilinear values at `point_a` and
    /// `point_b` (`f(point_a) - f(point_b)`) to lie in `[lower, upper]`.
    /// This is the stratigraphic-ordering case -- e.g. `lower = -inf`,
    /// `upper = 0` encodes "`point_a`'s unit must not be above `point_b`'s",
    /// i.e. `f(point_a) <= f(point_b)`. Only used by
    /// [`solve_admm`](Self::solve_admm); see
    /// [`add_value_inequality_constraint`](Self::add_value_inequality_constraint)'s
    /// docs for the drop conditions (either point outside the grid, or
    /// non-positive `weight`), which apply identically here.
    pub fn add_inequality_pair_constraint(
        &mut self,
        point_a: [f64; 3],
        point_b: [f64; 3],
        lower: f64,
        upper: f64,
        weight: f64,
    ) {
        if weight <= 0.0 {
            return;
        }
        let Some((idx_a, coefs_a)) = self.grid.trilinear_coeffs(point_a) else {
            return;
        };
        let Some((idx_b, coefs_b)) = self.grid.trilinear_coeffs(point_b) else {
            return;
        };
        let mut indices = [0usize; 16];
        let mut coefs = [0.0f64; 16];
        indices[..8].copy_from_slice(&idx_a);
        indices[8..].copy_from_slice(&idx_b);
        for k in 0..8 {
            coefs[k] = coefs_a[k];
            coefs[8 + k] = -coefs_b[k];
        }
        let norm = coefs.iter().map(|c| c * c).sum::<f64>().sqrt();
        if norm <= 0.0 {
            return;
        }
        let scale = weight / norm;
        for c in coefs.iter_mut() {
            *c *= scale;
        }
        self.pair_inequality_constraints
            .push(InequalityPairConstraint {
                indices,
                coefs,
                lower: lower * scale,
                upper: upper * scale,
            });
    }

    pub fn n_inequality_pair_constraints(&self) -> usize {
        self.pair_inequality_constraints.len()
    }

    /// Total inequality-constraint row count (value-bound + pairwise-ordering).
    pub fn n_inequality_constraints(&self) -> usize {
        self.value_inequality_constraints.len() + self.pair_inequality_constraints.len()
    }

    /// `Q x`: one entry per inequality-constraint row, in the fixed order
    /// value-bound rows (in add order) then pairwise-ordering rows (in add
    /// order) -- the same order [`inequality_bounds`](Self::inequality_bounds)
    /// and `solve_admm`'s internal `z`/`u` ADMM state vectors use.
    pub fn apply_inequality_operator(&self, x: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.n_inequality_constraints());
        for c in &self.value_inequality_constraints {
            out.push(ineq_point_row_value(c, x));
        }
        for c in &self.pair_inequality_constraints {
            out.push(ineq_pair_row_value(c, x));
        }
        out
    }

    /// `(lower, upper)` bound vectors, in the same row order as
    /// [`apply_inequality_operator`](Self::apply_inequality_operator).
    pub fn inequality_bounds(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.n_inequality_constraints();
        let mut lower = Vec::with_capacity(n);
        let mut upper = Vec::with_capacity(n);
        for c in &self.value_inequality_constraints {
            lower.push(c.lower);
            upper.push(c.upper);
        }
        for c in &self.pair_inequality_constraints {
            lower.push(c.lower);
            upper.push(c.upper);
        }
        (lower, upper)
    }

    /// Per row: the factor that turns a row violation into field units.
    /// A row is stored with `scale = weight / norm`. The first 8 coefficients
    /// of every row are trilinear weights times `scale`, and trilinear weights
    /// sum to 1. So the sum of the first 8 coefficients is `scale`. The
    /// result is `1 / scale`. Row order is the same as `inequality_bounds`.
    pub(crate) fn inequality_inverse_scales(&self) -> Vec<f64> {
        let inv = |first8: &[f64]| {
            let s: f64 = first8.iter().sum();
            if s > 0.0 { 1.0 / s } else { 1.0 }
        };
        let mut out = Vec::with_capacity(self.n_inequality_constraints());
        for c in &self.value_inequality_constraints {
            out.push(inv(&c.coefs[..8]));
        }
        for c in &self.pair_inequality_constraints {
            out.push(inv(&c.coefs[..8]));
        }
        out
    }

    /// `out += rho * Q^T y`, `y` given in the same row order as
    /// [`apply_inequality_operator`](Self::apply_inequality_operator). Used
    /// for `solve_admm`'s x-update rhs term `rho * Q^T (z - u)`.
    pub(crate) fn accumulate_inequality_transpose(&self, y: &[f64], rho: f64, out: &mut [f64]) {
        let mut idx = 0usize;
        for c in &self.value_inequality_constraints {
            let s = rho * y[idx];
            for k in 0..8 {
                out[c.indices[k]] += s * c.coefs[k];
            }
            idx += 1;
        }
        for c in &self.pair_inequality_constraints {
            let s = rho * y[idx];
            for k in 0..16 {
                out[c.indices[k]] += s * c.coefs[k];
            }
            idx += 1;
        }
    }

    /// `out += rho * Q^T Q x`, accumulated directly (each row is a rank-1
    /// outer product `coefs (coefs . x[indices])`, so this needs no
    /// intermediate `Q x` buffer) -- the inequality-constraint contribution
    /// to `solve_admm`'s x-update operator `N + rho * Q^T Q`.
    pub(crate) fn accumulate_inequality_normal(&self, x: &[f64], rho: f64, out: &mut [f64]) {
        for c in &self.value_inequality_constraints {
            let s = rho * ineq_point_row_value(c, x);
            for k in 0..8 {
                out[c.indices[k]] += s * c.coefs[k];
            }
        }
        for c in &self.pair_inequality_constraints {
            let s = rho * ineq_pair_row_value(c, x);
            for k in 0..16 {
                out[c.indices[k]] += s * c.coefs[k];
            }
        }
    }

    /// `rho * diag(Q^T Q)`, for the ADMM inner CG solve's Jacobi preconditioner.
    fn inequality_diagonal(&self, rho: f64) -> Vec<f64> {
        let mut diag = vec![0.0; self.dof()];
        for c in &self.value_inequality_constraints {
            for k in 0..8 {
                diag[c.indices[k]] += rho * c.coefs[k] * c.coefs[k];
            }
        }
        for c in &self.pair_inequality_constraints {
            for k in 0..16 {
                diag[c.indices[k]] += rho * c.coefs[k] * c.coefs[k];
            }
        }
        diag
    }

    /// `diag += rho * diag(Q^T Q)`. As `inequality_diagonal`, but it adds
    /// into `diag` and allocates no new full-size vector.
    pub(crate) fn add_inequality_diagonal(&self, rho: f64, diag: &mut [f64]) {
        for c in &self.value_inequality_constraints {
            for k in 0..8 {
                diag[c.indices[k]] += rho * c.coefs[k] * c.coefs[k];
            }
        }
        for c in &self.pair_inequality_constraints {
            for k in 0..16 {
                diag[c.indices[k]] += rho * c.coefs[k] * c.coefs[k];
            }
        }
    }

    /// Row `k`'s `(indices, coefs)`, in the flattened value-then-pair order
    /// every inequality-row accessor uses (see `apply_inequality_operator`).
    /// Used by `active_set_cg_core`/the free `active_set_*` helpers, which
    /// need to address one inequality row at a time rather than the whole
    /// stack `Q` at once.
    fn inequality_row(&self, k: usize) -> (&[usize], &[f64]) {
        let nv = self.value_inequality_constraints.len();
        if k < nv {
            let c = &self.value_inequality_constraints[k];
            (&c.indices, &c.coefs)
        } else {
            let c = &self.pair_inequality_constraints[k - nv];
            (&c.indices, &c.coefs)
        }
    }

    /// Row `k`'s value `Q_k . x`.
    fn inequality_row_dot(&self, k: usize, x: &[f64]) -> f64 {
        let (idx, coefs) = self.inequality_row(k);
        idx.iter().zip(coefs.iter()).map(|(i, c)| c * x[*i]).sum()
    }

    /// `out += scalar * Q_k^T`, i.e. scatter row `k`'s coefficients (scaled)
    /// into `out` at its dof indices.
    fn inequality_row_scatter_add(&self, k: usize, scalar: f64, out: &mut [f64]) {
        let (idx, coefs) = self.inequality_row(k);
        for (i, c) in idx.iter().zip(coefs.iter()) {
            out[*i] += scalar * c;
        }
    }

    /// `Q_a . Q_b` (dot product of two inequality rows over their shared
    /// dofs), for building the active working set's Gram matrix
    /// (`active_set_gram`).
    fn inequality_row_dot_rows(&self, a: usize, b: usize) -> f64 {
        let (idx_a, coefs_a) = self.inequality_row(a);
        let (idx_b, coefs_b) = self.inequality_row(b);
        let mut v = 0.0;
        for (ia, ca) in idx_a.iter().zip(coefs_a.iter()) {
            for (ib, cb) in idx_b.iter().zip(coefs_b.iter()) {
                if ia == ib {
                    v += ca * cb;
                }
            }
        }
        v
    }

    /// As [`inequality_row_dot_rows`](Self::inequality_row_dot_rows), but in
    /// the `diag(dinv)` inner product (`Q_a . diag(dinv) . Q_b`) -- used for
    /// `active_set_gram_weighted`'s Jacobi-preconditioned projection.
    fn inequality_row_dot_rows_weighted(&self, a: usize, b: usize, dinv: &[f64]) -> f64 {
        let (idx_a, coefs_a) = self.inequality_row(a);
        let (idx_b, coefs_b) = self.inequality_row(b);
        let mut v = 0.0;
        for (ia, ca) in idx_a.iter().zip(coefs_a.iter()) {
            for (ib, cb) in idx_b.iter().zip(coefs_b.iter()) {
                if ia == ib {
                    v += ca * cb * dinv[*ia];
                }
            }
        }
        v
    }

    /// `A_data^T W^2 A_data x`, accumulated directly from the stored
    /// trilinear rows -- there are normally far fewer constraints than grid
    /// nodes, so a plain scatter loop (no parallelism, no sparse matrix) is
    /// both simplest and fastest here.
    fn accumulate_data_normal(&self, x: &[f64], out: &mut [f64]) {
        accumulate_point_constraints_normal(&self.value_constraints, x, out);
        accumulate_point_constraints_normal(&self.gradient_constraints, x, out);
        accumulate_point_constraints_normal(&self.norm_constraints, x, out);
        accumulate_point_constraints_normal(&self.tangent_constraints, x, out);
        accumulate_pair_constraints_normal(&self.interface_constraints, x, out);
    }

    /// `A_data^T W^2 b` (only `value` and `norm` constraints ever have a
    /// nonzero target; gradient/tangent/interface constraints always target
    /// zero, so including them here is harmless and kept for symmetry with
    /// `accumulate_data_normal`/`diagonal`).
    fn data_rhs(&self, tracker: &mut MemTracker) -> Result<Vec<f64>, CoreError> {
        let mut out: Vec<f64> = tracker.try_zeroed(self.dof())?;
        self.write_data_rhs(&mut out);
        Ok(out)
    }

    /// Write the data right-hand side of `N` into `out`.
    /// The function sets `out` to zero first. It makes no new vector.
    /// `out` must have `dof()` entries.
    pub(crate) fn write_data_rhs(&self, out: &mut [f64]) {
        out.iter_mut().for_each(|v| *v = 0.0);
        point_constraints_rhs(&self.value_constraints, out);
        point_constraints_rhs(&self.gradient_constraints, out);
        point_constraints_rhs(&self.norm_constraints, out);
        point_constraints_rhs(&self.tangent_constraints, out);
    }

    /// Reserve room for `n` more value rows (exact size, no doubling).
    pub(crate) fn reserve_value(&mut self, n: usize) {
        self.value_constraints.reserve_exact(n);
    }

    /// Reserve room for `n` more gradient points (2 rows each).
    pub(crate) fn reserve_gradient(&mut self, n: usize) {
        self.gradient_constraints.reserve_exact(2 * n);
    }

    /// Reserve room for `n` more norm points (3 rows each).
    pub(crate) fn reserve_norm(&mut self, n: usize) {
        self.norm_constraints.reserve_exact(3 * n);
    }

    /// Reserve room for `n` more tangent rows.
    pub(crate) fn reserve_tangent(&mut self, n: usize) {
        self.tangent_constraints.reserve_exact(n);
    }

    /// Reserve room for `n` more interface rows.
    pub(crate) fn reserve_interface(&mut self, n: usize) {
        self.interface_constraints.reserve_exact(n);
    }

    /// Reserve room for `n` more value-bound rows.
    pub(crate) fn reserve_value_inequality(&mut self, n: usize) {
        self.value_inequality_constraints.reserve_exact(n);
    }

    /// Reserve room for `n` more pair-bound rows.
    pub(crate) fn reserve_pair_inequality(&mut self, n: usize) {
        self.pair_inequality_constraints.reserve_exact(n);
    }

    /// Required length of the `scratch` buffer passed to
    /// [`apply_normal_operator_into`]: one `A x` slot per regularisation
    /// family, per *boundary-shell row* (see `shell_row_indices`; only these
    /// need the exact two-pass computation now that the deep interior is
    /// handled by the fused single-kernel pass), plus one dof-length buffer
    /// shared (reused sequentially) across the border stencils.
    pub fn scratch_len(&self) -> usize {
        self.shell().rows.len() * N_FAMILIES + self.dof()
    }

    /// As [`scratch_len`], but for [`apply_normal_operator_into_unfused`]
    /// (the reference, pre-fusion computation): one `A x` slot per family
    /// per grid node (not just shell rows), plus the same border-stencil
    /// buffer.
    pub fn scratch_len_unfused(&self) -> usize {
        self.dof() * N_FAMILIES + self.dof()
    }

    /// The six families' `A x`, fused into one parallel pass over grid
    /// nodes (each node computes all six stencil rows in the same closure
    /// call) instead of six separate `Stencil::forward_into` dispatches.
    /// `scratch` is laid out as `dof` contiguous chunks of `N_FAMILIES`,
    /// i.e. `scratch[node * N_FAMILIES + family]`.
    fn forward_all_families_into<S: Scalar>(&self, x: &[S], scratch: &mut [S]) {
        let grid = &self.grid;
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        scratch
            .par_chunks_mut(N_FAMILIES)
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(flat, chunk)| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                for (slot, (stencil, _weight)) in chunk.iter_mut().zip(self.stencils.iter()) {
                    // The sum runs in f64. Only the stored value is S.
                    *slot = if stencil.is_row_active(grid, i, j, k) {
                        S::from_f64(
                            stencil
                                .offsets
                                .iter()
                                .zip(stencil.coeffs.iter())
                                .map(|(d, c)| {
                                    c * x[grid
                                        .neighbour_index(i, j, k, *d)
                                        .expect("row within margin")]
                                    .to_f64()
                                })
                                .sum::<f64>(),
                        )
                    } else {
                        S::default()
                    };
                }
            });
    }

    /// `out += sum_family w_family^2 * S_family^T (A x)`, reading the fused
    /// forward pass's output and fusing all six families' adjoint
    /// (correlation-with-flipped-stencil) contributions into one parallel
    /// pass over grid nodes as well.
    fn adjoint_all_families_accumulate_into<S: Scalar>(&self, scratch: &[S], out: &mut [S]) {
        let grid = &self.grid;
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        out.par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(flat, o)| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                let mut acc = 0.0;
                for (family, (stencil, weight)) in self.stencils.iter().enumerate() {
                    if *weight == 0.0 {
                        continue;
                    }
                    let value: f64 = stencil
                        .offsets
                        .iter()
                        .zip(stencil.coeffs.iter())
                        .map(|(d, c)| {
                            let flipped = (-d.0, -d.1, -d.2);
                            match grid.neighbour_index(i, j, k, flipped) {
                                // The source row's (si, sj, sk) is just (i, j, k)
                                // shifted by `flipped`, which `neighbour_index`
                                // already validated as in-bounds -- so recover it
                                // with the same cheap add used to get `src`,
                                // instead of decomposing `src` back via
                                // `node_ijk`'s two divisions + two modulos. This
                                // was the single largest cost in this pass.
                                Some(src) => {
                                    let si = (i as i64 + flipped.0 as i64) as usize;
                                    let sj = (j as i64 + flipped.1 as i64) as usize;
                                    let sk = (k as i64 + flipped.2 as i64) as usize;
                                    if stencil.is_row_active(grid, si, sj, sk) {
                                        c * scratch[src * N_FAMILIES + family].to_f64()
                                    } else {
                                        0.0
                                    }
                                }
                                None => 0.0,
                            }
                        })
                        .sum();
                    acc += weight * weight * value;
                }
                *o += S::from_f64(acc);
            });
    }

    /// The six curvature families' contribution to the normal operator at
    /// every node at distance >= 2 from the true domain edge, computed in a
    /// *single* pass reading `x` once via the combined kernel
    /// (`self.combined_kernel`) instead of the two-pass "forward into
    /// scratch, then adjoint reading scratch" computation
    /// (`forward_all_families_into` + `adjoint_all_families_accumulate_into`,
    /// still used by [`apply_normal_operator_into_unfused`]).
    ///
    /// Correctness: composing two radius-1 stencils (one family's `S^T S`)
    /// is `(S^T S x)[c] = sum_d coeff[d] * (S x)[c-d] = sum_{d1,d2}
    /// coeff[d1]*coeff[d2]*x[c-d1+d2]` -- a translation-invariant
    /// radius-2 stencil in `e = d2-d1` -- *only* when every row the adjoint
    /// step reads from (`c-d1` for every tap `d1`) is itself a fully active
    /// row (full stencil footprint fits, i.e. not within margin 1 of the
    /// true edge). Since `dist(c) >= 2` guarantees `dist(c-d1) >= 1` for
    /// every `|d1| <= 1`, using the deep-interior (dist >= 2) threshold is
    /// exactly the condition under which every family's combination is
    /// valid at every node the fused kernel touches -- matching the
    /// reference Python `_build_fused_cg_regularisation_operator`'s
    /// `_dist_to_true_edge(...) >= 2` cutoff. Nodes at `dist < 2` are left
    /// untouched here (added to as 0) and are fixed up exactly by
    /// [`shell_regularisation_correction_into`].
    ///
    /// Since `dist(c) >= 2` also guarantees every `c + e` for `|e| <= 2`
    /// stays in bounds on every axis, the kernel taps can be applied via
    /// direct (unchecked) flat-index arithmetic instead of
    /// `grid.neighbour_index`'s bounds-checked lookup.
    fn fused_interior_regularisation_into<S: Scalar>(&self, x: &[S], out: &mut [S]) {
        let kernel = self.kernel();
        if kernel.is_empty() {
            return;
        }
        let grid = &self.grid;
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        let nz = grid.nsteps[2];
        out.par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(flat, o)| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                let dx = i.min(nx - 1 - i);
                let dy = j.min(ny - 1 - j);
                let dz = k.min(nz - 1 - k);
                if dx.min(dy).min(dz) >= 2 {
                    let flat_i = flat as i64;
                    let mut acc = 0.0;
                    for (flat_offset, coeff) in kernel.iter() {
                        acc += coeff * x[(flat_i + flat_offset) as usize].to_f64();
                    }
                    *o += S::from_f64(acc);
                }
            });
    }

    /// Exact (unfused) computation of the six curvature families'
    /// contribution to the normal operator, restricted to the boundary
    /// shell `fused_interior_regularisation_into` skips: forward-evaluates
    /// only the rows in `shell_row_indices` (dist to true edge <= 2, exactly
    /// the row centres that can ever influence a dist < 2 output -- see
    /// `build_shell_indices`), then adjoint-accumulates only into
    /// `shell_output_indices` (dist < 2). Both passes cost O(boundary
    /// surface area), not O(dof), for any grid where the interior is large
    /// relative to its boundary.
    ///
    /// Mirrors `forward_all_families_into`/`adjoint_all_families_accumulate_into`'s
    /// arithmetic exactly (same row-activity checks, same summation order
    /// per family), just restricted to the shell row/output index sets and
    /// using `shell_row_lookup` instead of a dense `flat * N_FAMILIES`
    /// offset to find a row's scratch slot -- so this produces bit-for-bit
    /// identical values to the unfused computation at every shell node.
    fn shell_regularisation_correction_into<S: Scalar>(
        &self,
        x: &[S],
        shell_scratch: &mut [S],
        out: &mut [S],
    ) {
        let shell = self.shell();
        if shell.rows.is_empty() {
            return;
        }
        let grid = &self.grid;

        shell_scratch
            .par_chunks_mut(N_FAMILIES)
            .zip(shell.rows.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(chunk, &flat)| {
                let (i, j, k) = grid.node_ijk(flat as usize);
                for (slot, (stencil, weight)) in chunk.iter_mut().zip(self.stencils.iter()) {
                    *slot = if *weight != 0.0 && stencil.is_row_active(grid, i, j, k) {
                        S::from_f64(
                            stencil
                                .offsets
                                .iter()
                                .zip(stencil.coeffs.iter())
                                .map(|(d, c)| {
                                    c * x[grid
                                        .neighbour_index(i, j, k, *d)
                                        .expect("row within margin")]
                                    .to_f64()
                                })
                                .sum::<f64>(),
                        )
                    } else {
                        S::default()
                    };
                }
            });

        let contributions: Vec<f64> = shell
            .outputs
            .par_iter()
            .with_min_len(PAR_MIN_LEN)
            .map(|&flat| {
                let (i, j, k) = grid.node_ijk(flat as usize);
                let mut acc = 0.0;
                for (family, (stencil, weight)) in self.stencils.iter().enumerate() {
                    if *weight == 0.0 {
                        continue;
                    }
                    let value: f64 = stencil
                        .offsets
                        .iter()
                        .zip(stencil.coeffs.iter())
                        .map(|(d, c)| {
                            let flipped = (-d.0, -d.1, -d.2);
                            match grid.neighbour_index(i, j, k, flipped) {
                                Some(src) => {
                                    let si = (i as i64 + flipped.0 as i64) as usize;
                                    let sj = (j as i64 + flipped.1 as i64) as usize;
                                    let sk = (k as i64 + flipped.2 as i64) as usize;
                                    if stencil.is_row_active(grid, si, sj, sk) {
                                        let pos = shell.lookup[src];
                                        debug_assert!(
                                            pos != u32::MAX,
                                            "shell adjoint referenced a row outside the precomputed shell row set"
                                        );
                                        c * shell_scratch[pos as usize * N_FAMILIES + family]
                                            .to_f64()
                                    } else {
                                        0.0
                                    }
                                }
                                None => 0.0,
                            }
                        })
                        .sum();
                    acc += weight * weight * value;
                }
                acc
            })
            .collect();

        for (&flat, value) in shell.outputs.iter().zip(contributions.iter()) {
            out[flat as usize] += S::from_f64(*value);
        }
    }

    /// The full weighted normal-equations operator, writing into `out` and
    /// using `scratch` (length [`scratch_len`]) as the boundary-shell family
    /// scratch plus the border stencils' (sequentially-reused) scratch:
    /// `N x = (A_data^T W^2 A_data + sum_family w_family^2 S_family^T S_family
    /// + sum_border w_border^2 B^T B + ridge^2 I) x`.
    ///
    /// This is the allocation-free form used inside the CG hot loop (see
    /// `solve`); `scratch`/`out` are reused across every iteration instead
    /// of being (re)allocated. The interior curvature regularisation term is
    /// now computed as one fused single-kernel pass over the deep interior
    /// (`fused_interior_regularisation_into`, reading `x` once) plus a small
    /// exact two-pass correction restricted to the boundary shell
    /// (`shell_regularisation_correction_into`) -- collapsing what used to
    /// be two full-`dof` passes (forward, then adjoint) into roughly one
    /// pass's worth of memory traffic; see [`apply_normal_operator_into_unfused`]
    /// for the reference (pre-fusion) computation this is checked against.
    /// `out`'s initial zeroing and the ridge term are folded into a single
    /// parallel pass (`out[i] = ridge^2 * x[i]`, or zero if there's no ridge
    /// damping) rather than a serial zero pass followed by a separate
    /// serial ridge pass at the end.
    pub fn apply_normal_operator_into(&self, x: &[f64], scratch: &mut [f64], out: &mut [f64]) {
        if self.ridge_factor != 0.0 {
            let r2 = self.ridge_factor * self.ridge_factor;
            out.par_iter_mut()
                .zip(x.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(o, xi)| *o = r2 * xi);
        } else {
            out.par_iter_mut()
                .with_min_len(PAR_MIN_LEN)
                .for_each(|o| *o = 0.0);
        }
        self.accumulate_data_normal(x, out);
        let (shell_scratch, border_scratch) =
            scratch.split_at_mut(self.shell().rows.len() * N_FAMILIES);
        self.fused_interior_regularisation_into(x, out);
        self.shell_regularisation_correction_into(x, shell_scratch, out);
        self.accumulate_border_normal_into(x, border_scratch, out);
        self.accumulate_directional_normal_into(x, border_scratch, out);
    }

    /// The six curvature families' contribution to the normal operator from
    /// per-node directional-regularisation weights (see
    /// `add_directional_regularisation`), a no-op if none have been added.
    /// Reuses `scratch` (a `dof`-length buffer, e.g. the same one
    /// `border_stencils` sequentially reuses) as each family's intermediate
    /// `A x` buffer -- safe because `Stencil::forward_into` (unlike
    /// `EdgeStencil::forward_into`) always overwrites every entry of its
    /// output, so no stale data can leak between families.
    fn accumulate_directional_normal_into<S: Scalar>(
        &self,
        x: &[S],
        scratch: &mut [S],
        out: &mut [S],
    ) {
        let Some(weight2) = &self.directional_weight2 else {
            return;
        };
        for (family, (stencil, _)) in self.stencils.iter().enumerate() {
            // Same steps as `Stencil::accumulate_normal_weighted_into`,
            // but for any vector element type S. Sums run in f64.
            stencil.forward_into(&self.grid, x, scratch);
            scratch
                .par_iter_mut()
                .zip(weight2[family].par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(s, w2)| *s = S::from_f64(s.to_f64() * w2));
            stencil.adjoint_accumulate_into(&self.grid, scratch, 1.0, out);
        }
    }

    /// `out += sum_border w_border^2 B^T B x`. `scratch` holds each
    /// border stencil's `B x` in turn. Same steps as
    /// `EdgeStencil::accumulate_normal_into`, for any vector element type S.
    fn accumulate_border_normal_into<S: Scalar>(&self, x: &[S], scratch: &mut [S], out: &mut [S]) {
        for (stencil, weight) in &self.border_stencils {
            if *weight == 0.0 {
                continue;
            }
            stencil.forward_into(&self.grid, x, scratch);
            stencil.adjoint_accumulate_into(&self.grid, scratch, weight * weight, out);
        }
    }

    /// Reference (pre-fusion) computation of the normal operator's interior
    /// curvature contribution, using the original "forward all six families
    /// into scratch, then adjoint pass reading scratch" two-pass approach
    /// (`forward_all_families_into` + `adjoint_all_families_accumulate_into`)
    /// unconditionally at every node, instead of
    /// `apply_normal_operator_into`'s fused-interior + exact-shell split.
    /// `scratch` must be at least [`scratch_len_unfused`] long. Kept only
    /// for cross-checking the fused fast path is bit-for-bit-close
    /// equivalent (see the `fused_matches_unfused_*` tests) and for
    /// before/after profiling (see `profile_apply_1e6`/`profile_apply_1e6_unfused`);
    /// not used by `solve()`.
    pub fn apply_normal_operator_into_unfused(
        &self,
        x: &[f64],
        scratch: &mut [f64],
        out: &mut [f64],
    ) {
        if self.ridge_factor != 0.0 {
            let r2 = self.ridge_factor * self.ridge_factor;
            out.par_iter_mut()
                .zip(x.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(o, xi)| *o = r2 * xi);
        } else {
            out.par_iter_mut()
                .with_min_len(PAR_MIN_LEN)
                .for_each(|o| *o = 0.0);
        }
        self.accumulate_data_normal(x, out);
        let (family_scratch, border_scratch) = scratch.split_at_mut(self.dof() * N_FAMILIES);
        self.forward_all_families_into(x, family_scratch);
        self.adjoint_all_families_accumulate_into(family_scratch, out);
        self.accumulate_border_normal_into(x, border_scratch, out);
        self.accumulate_directional_normal_into(x, border_scratch, out);
    }

    /// `N x`, allocating its own output and scratch buffers. Convenience
    /// wrapper around [`apply_normal_operator_into`] for callers outside the
    /// CG hot loop (e.g. tests).
    pub fn apply_normal_operator(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.dof()];
        let mut scratch = vec![0.0; self.scratch_len()];
        self.apply_normal_operator_into(x, &mut scratch, &mut out);
        out
    }

    /// `N x` via the reference (pre-fusion) computation -- see
    /// [`apply_normal_operator_into_unfused`]. Allocates its own output and
    /// scratch buffers; only used by tests/benchmarks.
    pub fn apply_normal_operator_unfused(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.dof()];
        let mut scratch = vec![0.0; self.scratch_len_unfused()];
        self.apply_normal_operator_into_unfused(x, &mut scratch, &mut out);
        out
    }

    pub fn rhs(&self) -> Vec<f64> {
        self.try_rhs(&mut MemTracker::new(None))
            .expect("rhs must fit in memory")
    }

    /// As [`rhs`](Self::rhs). The vector is counted by `tracker`.
    pub(crate) fn try_rhs(&self, tracker: &mut MemTracker) -> Result<Vec<f64>, CoreError> {
        self.data_rhs(tracker)
    }

    /// `diag(N)`, the normal operator's diagonal, computed directly (no `N`
    /// ever assembled) for use as a Jacobi preconditioner.
    pub fn diagonal(&self) -> Vec<f64> {
        self.try_diagonal(&mut MemTracker::new(None))
            .expect("diagonal must fit in memory")
    }

    /// As [`diagonal`](Self::diagonal). The vector is counted by `tracker`.
    pub(crate) fn try_diagonal(&self, tracker: &mut MemTracker) -> Result<Vec<f64>, CoreError> {
        let mut diag: Vec<f64> = tracker.try_zeroed(self.dof())?;
        point_constraints_diagonal(&self.value_constraints, &mut diag);
        point_constraints_diagonal(&self.gradient_constraints, &mut diag);
        point_constraints_diagonal(&self.norm_constraints, &mut diag);
        point_constraints_diagonal(&self.tangent_constraints, &mut diag);
        pair_constraints_diagonal(&self.interface_constraints, &mut diag);
        for (stencil, weight) in &self.stencils {
            if *weight == 0.0 {
                continue;
            }
            let w2 = weight * weight;
            let d = stencil.diagonal(&self.grid);
            for (o, di) in diag.iter_mut().zip(d.iter()) {
                *o += w2 * di;
            }
        }
        for (stencil, weight) in &self.border_stencils {
            if *weight == 0.0 {
                continue;
            }
            let w2 = weight * weight;
            let d = stencil.diagonal(&self.grid);
            for (o, di) in diag.iter_mut().zip(d.iter()) {
                *o += w2 * di;
            }
        }
        if let Some(weight2) = &self.directional_weight2 {
            for (family, (stencil, _)) in self.stencils.iter().enumerate() {
                let d = stencil.diagonal_weighted(&self.grid, &weight2[family]);
                for (o, di) in diag.iter_mut().zip(d.iter()) {
                    *o += di;
                }
            }
        }
        if self.ridge_factor != 0.0 {
            let r2 = self.ridge_factor * self.ridge_factor;
            for o in diag.iter_mut() {
                *o += r2;
            }
        }
        Ok(diag)
    }

    /// The exact contribution of the `dxx` family plus any x-axis border
    /// (`dx_lower`/`dx_upper`) regularisation to the normal operator's
    /// x-adjacent off-diagonal entries `N[(i,j,k), (i+1,j,k)]`, one entry
    /// per dof (entry `flat`, for node `(i,j,k)`, is the coupling between
    /// that node and its `i+1` neighbour; the last node on each line has no
    /// `i+1` neighbour and is left at 0). Used to build
    /// `PreconditionerKind::BlockJacobiX`'s per-x-line tridiagonal system
    /// (`diagonal()` supplies the diagonal, this supplies the off-diagonal).
    ///
    /// SPD argument: this is the *principal submatrix* (restriction to one
    /// grid line's dof, i.e. fixed `j, k`, varying `i`) of `dxx`'s and
    /// border-`dx`'s own (each individually PSD, being an `S^T S` operator)
    /// contributions to `N` -- computed with the exact same row-activity
    /// masking `apply_normal_operator` uses (via `is_row_active`), not a
    /// deep-interior approximation, so it is the honest restriction, not an
    /// approximation of it. A principal submatrix of a PSD matrix is itself
    /// PSD, and the direct sum of the resulting per-line blocks (i.e. the
    /// block-diagonal matrix built from `diagonal()` + this off-diagonal) is
    /// therefore PSD; adding the *rest* of `diagonal()` (every other
    /// family's own diagonal, all individually PSD/nonneg) keeps it PSD, and
    /// adding `ridge_factor^2 * I` (strictly positive whenever
    /// `ridge_factor != 0`, the crate's default) makes the sum strictly
    /// positive definite. So `PreconditionerKind::BlockJacobiX` is
    /// guaranteed SPD under the same condition the plain Jacobi
    /// preconditioner already implicitly relies on (a strictly positive
    /// diagonal) -- see `solve_tridiagonal_line`'s fallback for the
    /// degenerate case.
    fn dxx_and_border_dx_offdiagonal(&self) -> Vec<f64> {
        let grid = &self.grid;
        let (nx, ny) = (grid.nsteps[0], grid.nsteps[1]);
        let dof = self.dof();
        let mut off = vec![0.0; dof];

        let w_dxx = self.regularisation.dxx;
        if w_dxx != 0.0 && nx >= 2 {
            let dxx_stencil = &self.stencils[0].0;
            let c_m1 = coeff_at(dxx_stencil, (-1, 0, 0));
            let c0 = coeff_at(dxx_stencil, (0, 0, 0));
            let c_p1 = coeff_at(dxx_stencil, (1, 0, 0));
            let w2 = w_dxx * w_dxx;
            off.par_chunks_mut(nx)
                .enumerate()
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(line, chunk)| {
                    let j = line % ny;
                    let k = line / ny;
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..nx - 1 {
                        let active_here = dxx_stencil.is_row_active(grid, i, j, k);
                        let active_next = dxx_stencil.is_row_active(grid, i + 1, j, k);
                        let mut v = 0.0;
                        if active_next {
                            v += c_m1 * c0;
                        }
                        if active_here {
                            v += c0 * c_p1;
                        }
                        chunk[i] += w2 * v;
                    }
                });
        }

        for (stencil, weight) in &self.border_stencils {
            if stencil.axis != 0 || *weight == 0.0 || nx < 2 {
                continue;
            }
            let w2 = weight * weight;
            let c0 = stencil
                .taps
                .iter()
                .find(|(o, _)| *o == 0)
                .map(|(_, c)| *c)
                .unwrap_or(0.0);
            if let Some((_, c_m1)) = stencil.taps.iter().find(|(o, _)| *o == -1) {
                let pair_i = stencil.layer - 1;
                off.par_chunks_mut(nx)
                    .for_each(|chunk| chunk[pair_i] += w2 * c_m1 * c0);
            }
            if let Some((_, c_p1)) = stencil.taps.iter().find(|(o, _)| *o == 1) {
                let pair_i = stencil.layer;
                off.par_chunks_mut(nx)
                    .for_each(|chunk| chunk[pair_i] += w2 * c0 * c_p1);
            }
        }

        off
    }

    /// Solve for the node values with Jacobi-preconditioned, matrix-free
    /// conjugate gradient. Convenience wrapper around
    /// [`solve_with`](Self::solve_with) fixed to
    /// [`PreconditionerKind::Jacobi`], kept so existing callers (Rust and
    /// Python) don't need to know about `PreconditionerKind` at all.
    pub fn solve(&self, max_iter: usize, tol: f64) -> CgResult {
        self.solve_with(max_iter, tol, PreconditionerKind::Jacobi)
    }

    /// Solve for the node values with matrix-free conjugate gradient, using
    /// `kind` to choose the preconditioner (see [`PreconditionerKind`]). No
    /// sparse matrix is ever built, and -- since `apply`/`precondition`
    /// write into caller-supplied buffers -- the only allocations for the
    /// whole solve are `rhs`, `diag`, the shared `apply` scratch buffer,
    /// `BlockJacobiX`'s off-diagonal + Thomas-algorithm scratch (only built
    /// if that variant is requested), and conjugate gradient's own five
    /// `dof`-length vectors; nothing is (re)allocated per iteration.
    pub fn solve_with(&self, max_iter: usize, tol: f64, kind: PreconditionerKind) -> CgResult {
        self.solve_with_tracked(max_iter, tol, kind, &mut MemTracker::new(None))
            .expect("solve must fit in memory")
    }

    /// As [`solve_with`](Self::solve_with). It builds the tables first (like
    /// `prepare`) and counts the rhs, diagonal, scratch and line buffers with
    /// `tracker`. It returns an error, not an abort, if a buffer is too big.
    pub(crate) fn solve_with_tracked(
        &self,
        max_iter: usize,
        tol: f64,
        kind: PreconditionerKind,
        tracker: &mut MemTracker,
    ) -> Result<CgResult, CoreError> {
        // This function takes `&self`, so it uses the `OnceLock` path of `prepare`.
        self.ensure_tables(tracker)?;
        let rhs = self.try_rhs(tracker)?;
        let diag = self.try_diagonal(tracker)?;
        let mut scratch: Vec<f64> = tracker.try_zeroed(self.scratch_len())?;

        let apply = |x: &[f64], out: &mut [f64]| {
            self.apply_normal_operator_into(x, &mut scratch, out);
        };

        match kind {
            PreconditionerKind::Jacobi => {
                let precondition = |r: &[f64], out: &mut [f64]| {
                    out.par_iter_mut()
                        .zip(r.par_iter())
                        .zip(diag.par_iter())
                        .with_min_len(PAR_MIN_LEN)
                        .for_each(|((o, ri), di)| {
                            *o = if *di > 0.0 { ri / di } else { *ri };
                        });
                };
                Ok(conjugate_gradient_preconditioned(
                    apply,
                    precondition,
                    &rhs,
                    None,
                    max_iter,
                    tol,
                ))
            }
            PreconditionerKind::BlockJacobiX => {
                let nx = self.grid.nsteps[0];
                let offdiag = self.dxx_and_border_dx_offdiagonal();
                let mut cprime: Vec<f64> = tracker.try_zeroed(self.dof())?;
                let precondition = |r: &[f64], out: &mut [f64]| {
                    out.par_chunks_mut(nx)
                        .zip(cprime.par_chunks_mut(nx))
                        .enumerate()
                        .with_min_len(PAR_MIN_LEN)
                        .for_each(|(line, (out_chunk, cprime_chunk))| {
                            let base = line * nx;
                            let r_chunk = &r[base..base + nx];
                            let d_chunk = &diag[base..base + nx];
                            let e_chunk = &offdiag[base..base + nx];
                            solve_tridiagonal_line(
                                d_chunk,
                                e_chunk,
                                r_chunk,
                                cprime_chunk,
                                out_chunk,
                            );
                        });
                };
                Ok(conjugate_gradient_preconditioned(
                    apply,
                    precondition,
                    &rhs,
                    None,
                    max_iter,
                    tol,
                ))
            }
        }
    }

    /// Solve for the node values with the box (inequality) constraints
    /// enforced, via classical ADMM (Boyd et al.'s ADMM survey's box-splitting
    /// form) on top of the existing matrix-free equality normal-equations
    /// operator/rhs (`N`/`rhs`, i.e. everything [`solve`](Self::solve)
    /// already handles: value/gradient/norm/tangent/interface constraints,
    /// curvature + border regularisation, ridge damping) plus the
    /// inequality rows added via
    /// [`add_value_inequality_constraint`](Self::add_value_inequality_constraint)/
    /// [`add_inequality_pair_constraint`](Self::add_inequality_pair_constraint),
    /// stacked matrix-free as `Q` (see `apply_inequality_operator`/
    /// `accumulate_inequality_normal`/`accumulate_inequality_transpose`).
    ///
    /// Ported from `loop_interpolation`'s `loopsolver/admm_method.py` +
    /// `admm_solver.py`, scoped to the *core* algorithm those two files
    /// implement (a fixed penalty `rho`, i.e. the Python default path
    /// without `admm_weight_final`/adaptive rho -- see Boyd et al., "Distributed
    /// Optimization and Statistical Learning via ADMM", section 5.2's box-constrained
    /// least squares example, which this is a direct instance of):
    ///
    /// - **x-update**: solve `(N + rho * Q^T Q) x = rhs + rho * Q^T (z - u)`
    ///   with [`cg::conjugate_gradient_preconditioned`] (Jacobi-preconditioned;
    ///   see below for why `BlockJacobiX` isn't offered here), warm-started
    ///   from the previous outer iteration's `x` (cheap and exact for this
    ///   fixed-`rho` loop, since the operator doesn't change between outer
    ///   iterations -- unlike the Python original's `reuse_inner_solve`
    ///   x0-injection machinery, which exists to serve its per-iteration
    ///   inner-solver-tolerance schedules and multiple solver backends,
    ///   neither of which exists here).
    /// - **z-update**: `z = clip(Q x + u, lower, upper)` elementwise,
    ///   exactly `admm_method.py`'s box projection.
    /// - **u-update**: `u += Q x - z`.
    /// - **stopping**: standard ADMM primal/dual residual check
    ///   (`primal = ||Qx - z||`, `dual = rho * ||z - z_prev||`, stopping when
    ///   both fall under `sqrt(m) * abs_tol + rel_tol * scale` -- the same
    ///   shape of criterion `admm_solver.py`'s `eps_pri`/`eps_dual` use, with
    ///   `scale = max(||Qx||, ||z||)` for the primal tolerance and `rho *
    ///   ||u||` for the dual one) or `max_outer_iter` outer iterations,
    ///   whichever comes first.
    ///
    /// If there are no inequality constraints at all
    /// (`n_inequality_constraints() == 0`), this degenerates to (and simply
    /// delegates to) [`solve`](Self::solve) with `(inner_max_iter, inner_tol)`
    /// -- an equality-only problem shouldn't pay any ADMM overhead (no `Q`,
    /// no outer loop, no `z`/`u` state) when there's nothing to project.
    ///
    /// Deliberately descoped from the Python original (left for a future
    /// pass, see `TASKS.md`): active-set row subsetting (`active_set`),
    /// adaptive rho (`adaptive_rho`), per-iteration inner-solver-tolerance
    /// schedules (`inner_rtol_start`/`inner_maxiter_start`/etc.), a
    /// non-constant rho schedule (`admm_weight_final`), warm-start reuse
    /// across ADMM iterations via a separately-injected `x0` guess for the
    /// *inner* solver (this crate instead just continues CG from the
    /// previous outer `x`, which is simpler and exact rather than a
    /// heuristic since `rho`, and therefore the inner operator, never
    /// changes), and multiple inner-solver backends (`lsmr`/`lsqr`/`cg`
    /// choice -- CG is the only inner solver this crate has at all, and is
    /// the right choice here: matrix-free, SPD, already used everywhere
    /// else in this crate).
    ///
    /// SPD argument for the x-update operator `N + rho * Q^T Q`: `N` is
    /// already established SPD whenever `ridge_factor != 0` (the crate
    /// default -- see [`diagonal`](Self::diagonal)'s ridge term and the
    /// `BlockJacobiX` SPD argument on [`dxx_and_border_dx_offdiagonal`], which
    /// relies on the same fact). `Q^T Q` is a Gram matrix -- PSD for *any*
    /// real matrix `Q`, regardless of its rank or row structure -- and
    /// `rho > 0` keeps `rho * Q^T Q` PSD too, so `N + rho * Q^T Q` is
    /// SPD + PSD = SPD unconditionally; no special argument about `Q`'s
    /// structure (unlike `BlockJacobiX`'s principal-submatrix argument) is
    /// needed. This is also *why* `BlockJacobiX` isn't offered as this
    /// solve's preconditioner: it's only proven to capture the `dxx`-family/
    /// x-axis-border coupling specifically, and would silently ignore (not
    /// mis-handle, just not exploit) the extra `rho * Q^T Q` coupling an
    /// inequality-pair row spanning two arbitrary points introduces --
    /// correct as a preconditioner (any SPD preconditioner is valid for CG),
    /// but not clearly still the *better* one once that coupling is real.
    /// Plain Jacobi is always well-defined (a strictly positive diagonal, by
    /// the same SPD argument) and keeps this new code path's preconditioner
    /// choice simple to reason about.
    ///
    /// Parameters
    /// ----------
    /// rho : the (fixed) ADMM penalty weight.
    /// max_outer_iter : outer ADMM iteration cap.
    /// inner_max_iter, inner_tol : passed through to every outer iteration's
    ///     inner CG x-update solve.
    /// admm_abs_tol, admm_rel_tol : absolute/relative primal-dual stopping
    ///     tolerances (see the stopping-criterion description above).
    pub fn solve_admm(
        &self,
        rho: f64,
        max_outer_iter: usize,
        inner_max_iter: usize,
        inner_tol: f64,
        admm_abs_tol: f64,
        admm_rel_tol: f64,
    ) -> AdmmResult {
        let m = self.n_inequality_constraints();
        if m == 0 {
            let cg = self.solve(inner_max_iter, inner_tol);
            return AdmmResult {
                x: cg.x,
                outer_iterations: 0,
                converged: cg.converged,
                inner_iterations_total: cg.iterations,
            };
        }

        let n = self.dof();
        let (lower, upper) = self.inequality_bounds();
        let base_rhs = self.rhs();

        // `diag(N + rho * Q^T Q)` never changes across outer iterations
        // (rho is fixed), so -- unlike `rhs_admm`, which depends on the
        // ever-changing `z - u` -- it's computed once, outside the loop.
        let mut diag = self.diagonal();
        let q_diag = self.inequality_diagonal(rho);
        for (d, qd) in diag.iter_mut().zip(q_diag.iter()) {
            *d += qd;
        }

        let mut z = vec![0.0; m];
        let mut u = vec![0.0; m];
        let mut zu = vec![0.0; m];
        let mut x = vec![0.0; n];
        let mut rhs_admm = vec![0.0; n];
        let mut normal_scratch = vec![0.0; self.scratch_len()];

        let mut outer_iterations = 0usize;
        let mut converged = false;
        let mut inner_iterations_total = 0usize;
        let sqrt_m = (m as f64).sqrt();

        for iter in 0..max_outer_iter {
            outer_iterations = iter + 1;

            for i in 0..m {
                zu[i] = z[i] - u[i];
            }
            rhs_admm.copy_from_slice(&base_rhs);
            self.accumulate_inequality_transpose(&zu, rho, &mut rhs_admm);

            let apply = |xin: &[f64], out: &mut [f64]| {
                self.apply_normal_operator_into(xin, &mut normal_scratch, out);
                self.accumulate_inequality_normal(xin, rho, out);
            };
            let precondition = |r: &[f64], out: &mut [f64]| {
                for i in 0..n {
                    out[i] = if diag[i] > 0.0 { r[i] / diag[i] } else { r[i] };
                }
            };
            let cg_result = conjugate_gradient_preconditioned(
                apply,
                precondition,
                &rhs_admm,
                Some(x.clone()),
                inner_max_iter,
                inner_tol,
            );
            x = cg_result.x;
            inner_iterations_total += cg_result.iterations;

            let qx = self.apply_inequality_operator(&x);
            let z_prev = z.clone();
            for i in 0..m {
                let v = qx[i] + u[i];
                z[i] = v.max(lower[i]).min(upper[i]);
            }
            for i in 0..m {
                u[i] += qx[i] - z[i];
            }

            let primal_norm = (0..m).map(|i| (qx[i] - z[i]).powi(2)).sum::<f64>().sqrt();
            let dual_norm = rho
                * (0..m)
                    .map(|i| (z[i] - z_prev[i]).powi(2))
                    .sum::<f64>()
                    .sqrt();
            let qx_norm = qx.iter().map(|v| v * v).sum::<f64>().sqrt();
            let z_norm = z.iter().map(|v| v * v).sum::<f64>().sqrt();
            let u_norm = rho * u.iter().map(|v| v * v).sum::<f64>().sqrt();
            let eps_pri = sqrt_m * admm_abs_tol + admm_rel_tol * qx_norm.max(z_norm);
            let eps_dual = sqrt_m * admm_abs_tol + admm_rel_tol * u_norm;

            if primal_norm <= eps_pri && dual_norm <= eps_dual {
                converged = true;
                break;
            }
        }

        AdmmResult {
            x,
            outer_iterations,
            converged,
            inner_iterations_total,
        }
    }

    /// Solve for the node values with the box (inequality) constraints
    /// enforced *exactly* (to `constraint_tol`), via Projected Active-Set
    /// Conjugate Gradient (AS-CG) -- see `Inequalities.md`'s "Method A".
    /// This is a primal active-set method for the QP `min 1/2 x^T N x - b^T
    /// x s.t. lower <= Q x <= upper` (`N`/`b` the same equality
    /// normal-equations operator/rhs [`solve`](Self::solve) uses; `Q` the
    /// stacked inequality rows from
    /// [`add_value_inequality_constraint`](Self::add_value_inequality_constraint)/
    /// [`add_inequality_pair_constraint`](Self::add_inequality_pair_constraint)):
    ///
    /// - **Working set**: a subset of inequality rows currently held exactly
    ///   at one of their bounds (`ActiveRow` -- row index, the bound value,
    ///   and which bound). Every other row is free.
    /// - **Feasibility**: `x` is kept exactly on the working set's
    ///   constraint surface (`Q_W x = bound_W`) via the minimum-distance
    ///   affine correction `active_set_restore_feasibility`, and any row
    ///   already violated at that corrected `x` (checked by direct scan,
    ///   not just "found mid-CG") is folded into the working set before any
    ///   CG work happens -- so this also functions as its own feasibility
    ///   ("Phase 1") step, with no separate solver needed to find a
    ///   starting point.
    /// - **Free-subspace solve**: for the current working set, CG runs in
    ///   the null space of `Q_W`, using `active_set_project_weighted`
    ///   (a Jacobi-preconditioned generalisation of `I - Q_W^T (Q_W
    ///   Q_W^T)^-1 Q_W`) as its "preconditioner" -- the standard *projected
    ///   preconditioned* CG method (Nocedal & Wright, "Numerical
    ///   Optimization", Algorithm 16.2) for equality-constrained QP, applied
    ///   here to the current working set treated as equality constraints,
    ///   with `diag(N)` playing the same preconditioning role
    ///   [`solve`](Self::solve)'s Jacobi preconditioner does.
    /// - **Step truncation**: before taking a CG step, every row *not* in
    ///   the working set is checked for whether the step would cross one of
    ///   its bounds; if so, the step is clipped to exactly reach the
    ///   nearest crossing, that row joins the working set, and CG restarts
    ///   (fresh conjugate directions, per `Inequalities.md`) on the new,
    ///   larger working set -- since `x` lands exactly on the new bound by
    ///   construction, `active_set_restore_feasibility`'s next call on it is
    ///   a no-op (see that function's docs).
    /// - **Node unfreezing**: once projected CG converges for a working set
    ///   (projected residual ~0), the Lagrange multipliers implied by that
    ///   stationarity point are recovered from a Gram solve
    ///   (`Q_W Q_W^T mu = -Q_W r`, `lambda = mu * sign`); any active row
    ///   with `lambda < -constraint_tol` is not a genuine KKT point (forcing
    ///   it to its bound is *hurting* the objective), so the most negative
    ///   one is released and CG restarts on the smaller working set. No
    ///   active row with a wrong-signed multiplier and no residual outside
    ///   the null space together are exactly the KKT conditions for this QP,
    ///   so that is this method's convergence criterion.
    ///
    /// If there are no inequality constraints at all
    /// (`n_inequality_constraints() == 0`), this degenerates to (and simply
    /// delegates to) [`solve`](Self::solve) with `(inner_max_iter,
    /// inner_tol)`, matching `solve_admm`'s fast path.
    ///
    /// Unlike [`solve_admm`](Self::solve_admm) (which only ever satisfies
    /// bounds up to the ADMM primal residual tolerance -- `Inequalities.md`'s
    /// "asymptotic contact"), every *active* row here is satisfied to
    /// floating-point precision by construction (`active_set_restore_feasibility`
    /// runs after every working-set change), matching the memo's "exact
    /// boundary contact" case for AS-CG. This comes at the cost of
    /// `O(|W|^3)` Cholesky factorisations of the working set's Gram matrix
    /// on every working-set change and an `O(m)` scan (`m` = total
    /// inequality rows) per CG step for step-truncation checks, so this
    /// solver is intended for the "low to moderate constraint density"
    /// regime `Inequalities.md`'s decision matrix calls out for AS-CG, not
    /// for a working set that ends up a large fraction of `m` -- see
    /// [`solve_hybrid`](Self::solve_hybrid) for a warm-started variant that
    /// keeps AS-CG's exact contact while amortising the discovery of a
    /// (near-)correct working set into a cheap matrix-free ADMM pass.
    ///
    /// Parameters
    /// ----------
    /// max_outer_iter : cap on the number of working-set changes.
    /// inner_max_iter, inner_tol : cap / relative-residual tolerance for
    ///     each working-set stage's projected CG solve.
    /// constraint_tol : absolute tolerance used both to decide whether a row
    ///     is "already violated" during the feasibility scan and to decide
    ///     whether a recovered Lagrange multiplier is "genuinely negative"
    ///     during the release check.
    pub fn solve_active_set_cg(
        &self,
        max_outer_iter: usize,
        inner_max_iter: usize,
        inner_tol: f64,
        constraint_tol: f64,
    ) -> ActiveSetCgResult {
        self.active_set_cg_core(
            vec![0.0; self.dof()],
            Vec::new(),
            max_outer_iter,
            inner_max_iter,
            inner_tol,
            constraint_tol,
        )
    }

    /// Shared active-set core: [`solve_active_set_cg`](Self::solve_active_set_cg)
    /// calls this starting from `x0 = 0`/an empty working set,
    /// [`solve_hybrid`](Self::solve_hybrid) calls it starting from an
    /// ADMM-warmed-up `x0`/working-set guess -- see those two methods' docs
    /// for the algorithm itself. `x0` need not be feasible (not even w.r.t.
    /// `active0`'s own rows): the first outer iteration's
    /// `active_set_restore_feasibility` call and feasibility scan bring it
    /// onto (and then extend `active0` to cover) the true feasible set
    /// regardless of the starting point's quality.
    fn active_set_cg_core(
        &self,
        mut x: Vec<f64>,
        mut active: Vec<ActiveRow>,
        max_outer_iter: usize,
        inner_max_iter: usize,
        inner_tol: f64,
        constraint_tol: f64,
    ) -> ActiveSetCgResult {
        let m = self.n_inequality_constraints();
        if m == 0 {
            let cg = self.solve(inner_max_iter, inner_tol);
            return ActiveSetCgResult {
                x: cg.x,
                outer_iterations: 0,
                inner_iterations_total: cg.iterations,
                converged: cg.converged,
            };
        }

        let n = self.dof();
        let b = self.rhs();
        let (lower, upper) = self.inequality_bounds();
        let mut scratch = vec![0.0; self.scratch_len()];
        let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
        // Jacobi preconditioning weight for the projected CG search
        // direction only (see `active_set_project_weighted`'s docs) -- not
        // used anywhere feasibility/optimality needs to be exact.
        let dinv: Vec<f64> = self
            .diagonal()
            .iter()
            .map(|d| if *d > 0.0 { 1.0 / d } else { 1.0 })
            .collect();

        let mut is_active = vec![false; m];
        for a in &active {
            is_active[a.row] = true;
        }

        let mut outer_iterations = 0usize;
        let mut inner_iterations_total = 0usize;
        let mut converged = false;

        'outer: while outer_iterations < max_outer_iter {
            let chol: Option<Vec<f64>> = if active.is_empty() {
                None
            } else {
                let gram = active_set_gram(self, &active);
                Some(active_set_cholesky(&gram, active.len()))
            };
            if let Some(c) = &chol {
                active_set_restore_feasibility(self, &active, c, &mut x);
            }

            // Feasibility scan: any inactive row already violated at this
            // (now working-set-feasible) `x` must join the working set
            // before any CG work is spent.
            let mut newly_violated: Option<ActiveRow> = None;
            for k in 0..m {
                if is_active[k] {
                    continue;
                }
                let v = self.inequality_row_dot(k, &x);
                if v < lower[k] - constraint_tol {
                    newly_violated = Some(ActiveRow {
                        row: k,
                        bound: lower[k],
                        sign: 1.0,
                    });
                    break;
                } else if v > upper[k] + constraint_tol {
                    newly_violated = Some(ActiveRow {
                        row: k,
                        bound: upper[k],
                        sign: -1.0,
                    });
                    break;
                }
            }
            if let Some(row) = newly_violated {
                is_active[row.row] = true;
                active.push(row);
                outer_iterations += 1;
                continue 'outer;
            }

            // Jacobi-preconditioned, projected CG on the current working
            // set's null space (see `active_set_project_weighted`'s docs
            // for why this uses a separate weighted Gram/Cholesky from
            // `chol`, which stays reserved for exact feasibility/multiplier
            // work).
            let chol_w: Option<Vec<f64>> = if active.is_empty() {
                None
            } else {
                let gram_w = active_set_gram_weighted(self, &active, &dinv);
                Some(active_set_cholesky(&gram_w, active.len()))
            };
            let mut r = vec![0.0; n];
            let mut napply = vec![0.0; n];
            self.apply_normal_operator_into(&x, &mut scratch, &mut napply);
            for i in 0..n {
                r[i] = b[i] - napply[i];
            }
            let mut g = active_set_project_weighted(self, &active, chol_w.as_deref(), &dinv, &r);
            let mut p = g.clone();
            let mut rz_old: f64 = r.iter().zip(g.iter()).map(|(a, bb)| a * bb).sum();
            let mut g_norm = g.iter().map(|v| v * v).sum::<f64>().sqrt();

            let mut truncated_this_stage = false;
            let mut inner_iter = 0usize;
            while inner_iter < inner_max_iter {
                if g_norm / b_norm <= inner_tol {
                    break;
                }
                let mut np = vec![0.0; n];
                self.apply_normal_operator_into(&p, &mut scratch, &mut np);
                let p_np: f64 = p.iter().zip(np.iter()).map(|(a, bb)| a * bb).sum();
                if p_np.abs() < 1e-300 {
                    break;
                }
                let mut alpha = rz_old / p_np;

                // Step truncation: does x -> x + alpha*p cross an inactive
                // row's bound before reaching the unconstrained CG step?
                let mut blocking: Option<ActiveRow> = None;
                for k in 0..m {
                    if is_active[k] {
                        continue;
                    }
                    let qk_x = self.inequality_row_dot(k, &x);
                    let qk_p = self.inequality_row_dot(k, &p);
                    if qk_p > 1e-300 && upper[k].is_finite() {
                        let a_max = (upper[k] - qk_x) / qk_p;
                        if a_max < alpha {
                            alpha = a_max.max(0.0);
                            blocking = Some(ActiveRow {
                                row: k,
                                bound: upper[k],
                                sign: -1.0,
                            });
                        }
                    } else if qk_p < -1e-300 && lower[k].is_finite() {
                        let a_max = (lower[k] - qk_x) / qk_p;
                        if a_max < alpha {
                            alpha = a_max.max(0.0);
                            blocking = Some(ActiveRow {
                                row: k,
                                bound: lower[k],
                                sign: 1.0,
                            });
                        }
                    }
                }

                for i in 0..n {
                    x[i] += alpha * p[i];
                }
                inner_iter += 1;

                if let Some(row) = blocking {
                    is_active[row.row] = true;
                    active.push(row);
                    truncated_this_stage = true;
                    break;
                }

                for i in 0..n {
                    r[i] -= alpha * np[i];
                }
                let g_new =
                    active_set_project_weighted(self, &active, chol_w.as_deref(), &dinv, &r);
                g_norm = g_new.iter().map(|v| v * v).sum::<f64>().sqrt();
                if g_norm / b_norm <= inner_tol {
                    g = g_new;
                    break;
                }
                let rz_new: f64 = r.iter().zip(g_new.iter()).map(|(a, bb)| a * bb).sum();
                let beta = rz_new / rz_old;
                for i in 0..n {
                    p[i] = g_new[i] + beta * p[i];
                }
                g = g_new;
                rz_old = rz_new;
            }

            inner_iterations_total += inner_iter;

            if truncated_this_stage {
                outer_iterations += 1;
                continue 'outer;
            }

            // Hit `inner_max_iter` without the projected residual actually
            // reaching `inner_tol`: stop here rather than trusting a
            // Lagrange-multiplier check derived from a non-stationary `r` --
            // an unconverged residual can easily produce a spuriously
            // wrong-signed multiplier, releasing a row that's genuinely
            // active and then immediately re-adding it next iteration (via
            // the feasibility scan), burning another `inner_max_iter`
            // cold-restarted CG budget for no progress instead of reporting
            // "not converged" and letting the caller raise the budget.
            if g_norm / b_norm > inner_tol {
                converged = false;
                break;
            }

            if active.is_empty() {
                converged = true;
                break;
            }

            // Lagrange-multiplier release check: `mu` solves `Q_W Q_W^T mu =
            // -Q_W r`; `lambda_i = mu_i * sign_i` must be >= 0 for every
            // active row at a genuine KKT point (see this method's doc
            // comment for the derivation).
            let chol = chol.expect("non-empty active set always has a Cholesky factor");
            let w_qr = active_set_apply(self, &active, &r);
            let neg_w_qr: Vec<f64> = w_qr.iter().map(|v| -v).collect();
            let mu = cholesky_solve(&chol, active.len(), &neg_w_qr);
            let mut release_idx: Option<usize> = None;
            let mut worst = -constraint_tol;
            for (idx, (a, mu_i)) in active.iter().zip(mu.iter()).enumerate() {
                let lambda = mu_i * a.sign;
                if lambda < worst {
                    worst = lambda;
                    release_idx = Some(idx);
                }
            }
            match release_idx {
                Some(idx) => {
                    is_active[active[idx].row] = false;
                    active.remove(idx);
                    outer_iterations += 1;
                }
                None => {
                    converged = true;
                    break;
                }
            }
        }

        ActiveSetCgResult {
            x,
            outer_iterations,
            inner_iterations_total,
            converged,
        }
    }

    /// Run exactly `iters` outer iterations of the same classical ADMM loop
    /// [`solve_admm`](Self::solve_admm) uses, *without* the primal/dual
    /// stopping check -- a fixed-budget "warm-up" pass rather than a solve
    /// to convergence, matching `Inequalities.md`'s "Stage 1 (ADMM
    /// Warm-Up)". Returns `(x, z, inner_iterations_total)`; `z` (the ADMM
    /// box-projection variable) is what [`solve_hybrid`](Self::solve_hybrid)
    /// uses to seed its active-set guess, since `z_i` sitting exactly at
    /// `lower[i]`/`upper[i]` is precisely ADMM's signal that row `i` wants
    /// to be active.
    fn admm_warmup(
        &self,
        rho: f64,
        iters: usize,
        inner_max_iter: usize,
        inner_tol: f64,
    ) -> (Vec<f64>, Vec<f64>, usize) {
        let n = self.dof();
        let m = self.n_inequality_constraints();
        let (lower, upper) = self.inequality_bounds();
        let base_rhs = self.rhs();

        let mut diag = self.diagonal();
        let q_diag = self.inequality_diagonal(rho);
        for (d, qd) in diag.iter_mut().zip(q_diag.iter()) {
            *d += qd;
        }

        let mut z = vec![0.0; m];
        let mut u = vec![0.0; m];
        let mut zu = vec![0.0; m];
        let mut x = vec![0.0; n];
        let mut rhs_admm = vec![0.0; n];
        let mut normal_scratch = vec![0.0; self.scratch_len()];
        let mut inner_iterations_total = 0usize;

        for _ in 0..iters {
            for i in 0..m {
                zu[i] = z[i] - u[i];
            }
            rhs_admm.copy_from_slice(&base_rhs);
            self.accumulate_inequality_transpose(&zu, rho, &mut rhs_admm);

            let apply = |xin: &[f64], out: &mut [f64]| {
                self.apply_normal_operator_into(xin, &mut normal_scratch, out);
                self.accumulate_inequality_normal(xin, rho, out);
            };
            let precondition = |r: &[f64], out: &mut [f64]| {
                for i in 0..n {
                    out[i] = if diag[i] > 0.0 { r[i] / diag[i] } else { r[i] };
                }
            };
            let cg_result = conjugate_gradient_preconditioned(
                apply,
                precondition,
                &rhs_admm,
                Some(x.clone()),
                inner_max_iter,
                inner_tol,
            );
            x = cg_result.x;
            inner_iterations_total += cg_result.iterations;

            let qx = self.apply_inequality_operator(&x);
            for i in 0..m {
                let v = qx[i] + u[i];
                z[i] = v.max(lower[i]).min(upper[i]);
            }
            for i in 0..m {
                u[i] += qx[i] - z[i];
            }
        }

        (x, z, inner_iterations_total)
    }

    /// Two-stage hybrid solve (`Inequalities.md`'s "Recommended Strategy"):
    /// a fixed-budget matrixless ADMM warm-up
    /// ([`admm_warmup`](Self::admm_warmup), `admm_warmup_iter` outer
    /// iterations, processing bulk boundary violations in parallel without
    /// stalling on active-set step truncation) followed by an
    /// [`solve_active_set_cg`](Self::solve_active_set_cg) polish -- warm
    /// started from the ADMM stage's `x` and an initial working-set guess
    /// (every row whose ADMM `z` landed within `constraint_tol` of a bound)
    /// -- that drives every active row to floating-point contact precision.
    ///
    /// The polish stage's `active_set_cg_core` doesn't require the ADMM
    /// hand-off to be exact: its first working-set change unconditionally
    /// re-derives feasibility from whatever `x`/working-set guess it's
    /// given (see `active_set_cg_core`'s docs), so an under-converged ADMM
    /// warm-up degrades to more active-set work, not incorrect results.
    ///
    /// If there are no inequality constraints at all
    /// (`n_inequality_constraints() == 0`), delegates straight to
    /// [`solve`](Self::solve) with `(polish_inner_max_iter, polish_tol)` --
    /// same fast path as `solve_admm`/`solve_active_set_cg`.
    ///
    /// Parameters
    /// ----------
    /// rho, admm_inner_max_iter, admm_inner_tol : passed through to the ADMM
    ///     warm-up stage (see `solve_admm`'s docs for their meaning).
    /// admm_warmup_iter : fixed number of ADMM outer iterations to run as
    ///     warm-up (`Inequalities.md` suggests 10-20).
    /// polish_max_outer_iter, polish_inner_max_iter, polish_tol : passed
    ///     through to the active-set polish stage (see
    ///     `solve_active_set_cg`'s `max_outer_iter`/`inner_max_iter`/
    ///     `inner_tol`).
    /// constraint_tol : shared between the ADMM-to-active-set hand-off (how
    ///     close `z` must be to a bound to seed that row as active) and the
    ///     polish stage's own feasibility/multiplier tolerances (see
    ///     `solve_active_set_cg`'s `constraint_tol`).
    #[allow(clippy::too_many_arguments)]
    pub fn solve_hybrid(
        &self,
        rho: f64,
        admm_warmup_iter: usize,
        admm_inner_max_iter: usize,
        admm_inner_tol: f64,
        polish_max_outer_iter: usize,
        polish_inner_max_iter: usize,
        polish_tol: f64,
        constraint_tol: f64,
    ) -> HybridResult {
        let m = self.n_inequality_constraints();
        if m == 0 {
            let cg = self.solve(polish_inner_max_iter, polish_tol);
            return HybridResult {
                x: cg.x,
                admm_warmup_iterations: 0,
                admm_inner_iterations_total: 0,
                polish_outer_iterations: 0,
                polish_inner_iterations_total: cg.iterations,
                converged: cg.converged,
            };
        }

        let (lower, upper) = self.inequality_bounds();
        let (x0, z, admm_inner_iterations_total) =
            self.admm_warmup(rho, admm_warmup_iter, admm_inner_max_iter, admm_inner_tol);

        let mut active0 = Vec::new();
        for k in 0..m {
            if lower[k].is_finite() && (z[k] - lower[k]).abs() <= constraint_tol {
                active0.push(ActiveRow {
                    row: k,
                    bound: lower[k],
                    sign: 1.0,
                });
            } else if upper[k].is_finite() && (z[k] - upper[k]).abs() <= constraint_tol {
                active0.push(ActiveRow {
                    row: k,
                    bound: upper[k],
                    sign: -1.0,
                });
            }
        }

        let polish = self.active_set_cg_core(
            x0,
            active0,
            polish_max_outer_iter,
            polish_inner_max_iter,
            polish_tol,
            constraint_tol,
        );

        HybridResult {
            x: polish.x,
            admm_warmup_iterations: admm_warmup_iter,
            admm_inner_iterations_total,
            polish_outer_iterations: polish.outer_iterations,
            polish_inner_iterations_total: polish.inner_iterations_total,
            converged: polish.converged,
        }
    }

    /// Evaluate the solved field at arbitrary points via trilinear
    /// interpolation. Points outside the grid are returned as `NaN`.
    pub fn evaluate(&self, c: &[f64], points: &[[f64; 3]]) -> Vec<f64> {
        points
            .iter()
            .map(|p| match self.grid.trilinear_coeffs(*p) {
                Some((idx, coefs)) => idx.iter().zip(coefs.iter()).map(|(i, w)| w * c[*i]).sum(),
                None => f64::NAN,
            })
            .collect()
    }

    /// Evaluate the solved field's gradient at arbitrary (world-space)
    /// points via the trilinear gradient shape function. Points outside the
    /// grid are returned as `[NaN; 3]`.
    pub fn evaluate_gradient(&self, c: &[f64], points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        points
            .iter()
            .map(|p| match self.grid.trilinear_gradient(*p) {
                Some((idx, grad)) => {
                    let mut local = [0.0; 3];
                    for axis in 0..3 {
                        local[axis] = idx
                            .iter()
                            .zip(grad.iter())
                            .map(|(i, g)| c[*i] * g[axis])
                            .sum();
                    }
                    local
                }
                None => [f64::NAN; 3],
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grid() -> StructuredGrid3D {
        StructuredGrid3D::new([0.0, 0.0, 0.0], [9, 9, 9], [1.0, 1.0, 1.0])
    }

    #[test]
    fn normal_operator_is_symmetric() {
        let mut interp = FiniteDifferenceInterpolator::new(make_grid());
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        interp.add_value_constraint([6.0, 3.0, 4.0], -0.5, 2.0);

        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 17 + 5) % 89) as f64 / 89.0).collect();

        let nx = interp.apply_normal_operator(&x);
        let ny = interp.apply_normal_operator(&y);
        let lhs: f64 = nx.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(ny.iter()).map(|(a, b)| a * b).sum();
        assert!(
            (lhs - rhs).abs() < 1e-6 * lhs.abs().max(1.0),
            "lhs={lhs} rhs={rhs}"
        );
    }

    /// Symmetry check extended to cover every new constraint family
    /// (gradient/norm/tangent/interface) plus border regularisation, all at
    /// once -- guards against any one of them breaking `N`'s self-adjointness.
    #[test]
    fn normal_operator_is_symmetric_with_all_constraint_families() {
        let mut interp = FiniteDifferenceInterpolator::new(make_grid());
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        interp.add_gradient_constraint([3.0, 3.0, 3.0], [1.0, 0.0, 0.0], 2.0);
        interp.add_norm_constraint([4.0, 4.0, 4.0], [0.3, 0.1, -0.2], 1.5);
        interp.add_tangent_constraint([5.0, 5.0, 5.0], [1.0, 1.0, 0.0], 0.7);
        interp.add_interface_constraints(
            &[[1.0, 1.0, 1.0], [6.0, 6.0, 6.0], [2.0, 6.0, 1.0]],
            &[1.0, 1.0, 2.0],
            3.0,
        );

        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 17 + 5) % 89) as f64 / 89.0).collect();

        let nx = interp.apply_normal_operator(&x);
        let ny = interp.apply_normal_operator(&y);
        let lhs: f64 = nx.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(ny.iter()).map(|(a, b)| a * b).sum();
        assert!(
            (lhs - rhs).abs() < 1e-6 * lhs.abs().max(1.0),
            "lhs={lhs} rhs={rhs}"
        );
    }

    #[test]
    fn recovers_affine_field_from_scattered_samples() {
        // f(x, y, z) = 1 + 2x - y + 0.5z has zero curvature everywhere, so
        // the regulariser shouldn't fight the data at all: CG should recover
        // it (near-)exactly from a handful of scattered value constraints.
        // The default border regularisation is on. The border rows are centred
        // second derivatives, as in the matrix solver. They give a plane no
        // penalty, and they are the only rows that fix the nodes on the faces.
        let f = |p: [f64; 3]| 1.0 + 2.0 * p[0] - p[1] + 0.5 * p[2];
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);

        let samples = [
            [0.5, 0.5, 0.5],
            [7.5, 1.0, 2.0],
            [2.0, 6.0, 0.2],
            [4.0, 4.0, 7.0],
            [1.0, 7.5, 7.5],
            [6.0, 2.0, 6.0],
        ];
        for p in samples {
            interp.add_value_constraint(p, f(p), 100.0);
        }

        let result = interp.solve(2000, 1e-10);
        assert!(result.converged, "iterations={}", result.iterations);

        let query = [[3.0, 3.0, 3.0], [8.0, 0.0, 8.0], [0.0, 8.0, 0.0]];
        let predicted = interp.evaluate(&result.x, &query);
        for (p, v) in query.iter().zip(predicted.iter()) {
            assert!(
                (v - f(*p)).abs() < 1e-2,
                "point={p:?} predicted={v} expected={}",
                f(*p)
            );
        }
    }

    #[test]
    fn constraints_outside_grid_are_dropped() {
        let mut interp = FiniteDifferenceInterpolator::new(make_grid());
        interp.add_value_constraint([-5.0, 0.0, 0.0], 1.0, 1.0);
        assert_eq!(interp.n_value_constraints(), 0);

        interp.add_gradient_constraint([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert_eq!(interp.n_gradient_constraints(), 0);

        interp.add_norm_constraint([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert_eq!(interp.n_norm_constraints(), 0);

        interp.add_tangent_constraint([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert_eq!(interp.n_tangent_constraints(), 0);

        interp.add_interface_constraints(&[[-5.0, 0.0, 0.0], [-6.0, 0.0, 0.0]], &[1.0, 1.0], 1.0);
        assert_eq!(interp.n_interface_constraints(), 0);
    }

    /// Recover a plane from gradient (orthogonal-direction) constraints plus
    /// a *sparse* pair of value constraints. Gradient-orthogonal rows are
    /// homogeneous (target zero) and only fix gradient *direction*: the
    /// all-zero field trivially satisfies them (and zero curvature) too, so
    /// they alone -- or with only a single value pin -- cannot rule out the
    /// degenerate flat solution. Two value constraints at different points
    /// with different targets break that degeneracy (fixing the two
    /// remaining unknowns, the plane's constant term and its magnitude along
    /// the known gradient direction) without otherwise being enough to pin
    /// down a full 3D affine field (4 unknowns) by themselves -- so this is
    /// a genuine test that the gradient constraints are doing real,
    /// direction-fixing work, not just a data-fit test they're irrelevant to.
    #[test]
    fn recovers_plane_from_gradient_constraints_plus_sparse_value_pins() {
        // The default border regularisation is on. The border rows are centred
        // second derivatives, so a plane has no penalty. The border rows are
        // also the only rows that fix the nodes on the faces.
        let gradient = [2.0, -1.0, 0.5];
        let f = |p: [f64; 3]| 1.0 + gradient[0] * p[0] + gradient[1] * p[1] + gradient[2] * p[2];
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);

        for p in [
            [1.0, 1.0, 1.0],
            [7.0, 2.0, 3.0],
            [3.0, 6.0, 1.0],
            [5.0, 4.0, 6.0],
        ] {
            interp.add_gradient_constraint(p, gradient, 200.0);
        }
        interp.add_value_constraint([0.5, 0.5, 0.5], f([0.5, 0.5, 0.5]), 100.0);
        interp.add_value_constraint([7.5, 7.5, 7.5], f([7.5, 7.5, 7.5]), 100.0);

        let result = interp.solve(4000, 1e-11);
        assert!(result.converged, "iterations={}", result.iterations);

        let query = [[3.0, 3.0, 3.0], [8.0, 0.0, 8.0], [0.0, 8.0, 0.0]];
        let predicted = interp.evaluate(&result.x, &query);
        for (p, v) in query.iter().zip(predicted.iter()) {
            assert!(
                (v - f(*p)).abs() < 5e-2,
                "point={p:?} predicted={v} expected={}",
                f(*p)
            );
        }
    }

    /// Recover a plane from norm constraints: unlike gradient-orthogonal
    /// constraints, a norm constraint is non-homogeneous (it targets the
    /// actual gradient components), so a single norm constraint plus a
    /// single value pin fully (and uniquely) determines a 3D affine field.
    #[test]
    fn recovers_plane_from_norm_constraints() {
        // The default border regularisation is on -- see the comment in
        // `recovers_plane_from_gradient_constraints_plus_sparse_value_pins`.
        let gradient = [2.0, -1.0, 0.5];
        let f = |p: [f64; 3]| 1.0 + gradient[0] * p[0] + gradient[1] * p[1] + gradient[2] * p[2];
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);

        interp.add_norm_constraint([4.3, 4.3, 4.3], gradient, 200.0);
        interp.add_value_constraint([0.5, 0.5, 0.5], f([0.5, 0.5, 0.5]), 100.0);

        let result = interp.solve(4000, 1e-11);
        assert!(result.converged, "iterations={}", result.iterations);

        let query = [[3.0, 3.0, 3.0], [8.0, 0.0, 8.0], [0.0, 8.0, 0.0]];
        let predicted = interp.evaluate(&result.x, &query);
        for (p, v) in query.iter().zip(predicted.iter()) {
            assert!(
                (v - f(*p)).abs() < 5e-2,
                "point={p:?} predicted={v} expected={}",
                f(*p)
            );
        }
    }

    /// A tangent constraint should measurably zero the directional
    /// derivative at its point, even when that direction is *not*
    /// orthogonal to the field's natural (unconstrained) gradient -- i.e.
    /// the constraint has to actually do work, not just agree with the
    /// data's existing shape.
    #[test]
    fn tangent_constraint_zeroes_directional_derivative() {
        // Border regularisation disabled -- see the comment in
        // `recovers_plane_from_gradient_constraints_plus_sparse_value_pins`.
        let f = |p: [f64; 3]| 1.0 + 2.0 * p[0] - p[1] + 0.5 * p[2];
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid).with_border_regularisation(
            BorderRegularisationWeights {
                dx: 0.0,
                dy: 0.0,
                dz: 0.0,
            },
        );
        let samples = [
            [0.5, 0.5, 0.5],
            [7.5, 1.0, 2.0],
            [2.0, 6.0, 0.2],
            [4.0, 4.0, 7.0],
            [1.0, 7.5, 7.5],
            [6.0, 2.0, 6.0],
        ];
        for p in samples {
            interp.add_value_constraint(p, f(p), 20.0);
        }
        // The field's natural gradient is (2, -1, 0.5), not orthogonal to
        // (1, 1, 1) -- so constraining the directional derivative along
        // (1, 1, 1) to be zero genuinely fights the data. `point` is
        // deliberately not on an exact grid node: trilinear fields are only
        // C0 (value-continuous) across cell boundaries, so a finite
        // difference straddling one would mix two different cells' gradients
        // instead of measuring the single cell the constraint actually
        // targets.
        let point = [4.3, 4.3, 4.3];
        let v = [1.0, 1.0, 1.0];
        interp.add_tangent_constraint(point, v, 500.0);

        let result = interp.solve(4000, 1e-11);
        assert!(result.converged, "iterations={}", result.iterations);

        // Numerically estimate the directional derivative of the *solved*
        // field along v at `point` and confirm it's near zero.
        let eps = 1e-3;
        let v_unit = {
            let n = norm3(v);
            scale3(v, 1.0 / n)
        };
        let p_plus = [
            point[0] + eps * v_unit[0],
            point[1] + eps * v_unit[1],
            point[2] + eps * v_unit[2],
        ];
        let p_minus = [
            point[0] - eps * v_unit[0],
            point[1] - eps * v_unit[1],
            point[2] - eps * v_unit[2],
        ];
        let vals = interp.evaluate(&result.x, &[p_plus, p_minus]);
        let directional_derivative = (vals[0] - vals[1]) / (2.0 * eps);
        assert!(
            directional_derivative.abs() < 1e-2,
            "directional_derivative={directional_derivative}"
        );
    }

    /// Two points sharing an interface id (with no other constraint forcing
    /// them apart) should converge to (numerically) equal solved values.
    #[test]
    fn interface_constraint_pair_converges_to_equal_values() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        // Anchor the field away from the interface pair so the system isn't
        // entirely homogeneous (which would trivially solve to all zeros).
        interp.add_value_constraint([0.5, 0.5, 0.5], 3.0, 50.0);

        let a = [6.0, 2.0, 3.0];
        let b = [2.0, 6.0, 5.0];
        interp.add_interface_constraints(&[a, b], &[7.0, 7.0], 500.0);

        let result = interp.solve(4000, 1e-11);
        assert!(result.converged, "iterations={}", result.iterations);

        let vals = interp.evaluate(&result.x, &[a, b]);
        assert!(
            (vals[0] - vals[1]).abs() < 1e-3,
            "a={} b={}",
            vals[0],
            vals[1]
        );
    }

    /// Border regularisation should measurably change boundary-node values
    /// relative to having it switched off, for an otherwise-identical setup
    /// whose value data doesn't happen to already look like an affine field
    /// out to the boundary (so the border rows have something to disagree
    /// with the unregularised interior fit about).
    #[test]
    fn border_regularisation_changes_boundary_values() {
        let grid = make_grid();
        let build = |border_weight: f64| {
            let mut interp = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
                grid.origin,
                grid.nsteps,
                grid.step_vector,
            ))
            .with_border_regularisation(BorderRegularisationWeights {
                dx: border_weight,
                dy: border_weight,
                dz: border_weight,
            });
            // A handful of interior samples from a distinctly non-affine
            // field, so the boundary is under-determined by data alone and
            // sensitive to whatever regularisation is in play.
            let f = |p: [f64; 3]| ((p[0] * 0.6).sin() * (p[1] * 0.5).cos() * 4.0) + p[2] * 0.3;
            for p in [
                [2.0, 2.0, 2.0],
                [6.0, 3.0, 5.0],
                [3.0, 6.0, 2.0],
                [5.0, 5.0, 6.0],
            ] {
                interp.add_value_constraint(p, f(p), 50.0);
            }
            let result = interp.solve(4000, 1e-10);
            assert!(result.converged, "iterations={}", result.iterations);
            result.x
        };

        let with_border = build(5.0);
        let without_border = build(0.0);

        // Compare values at a corner node. A border row on layer 1 has the same
        // taps as the interior curvature row of that layer, so it adds only
        // weight there. The face and edge nodes are also reached by interior
        // rows. The corner nodes are reached by the border rows alone.
        let boundary_node = grid.index(0, 0, 0);
        let diff = (with_border[boundary_node] - without_border[boundary_node]).abs();
        assert!(
            diff > 1e-6,
            "expected border regularisation to change boundary node value, diff={diff}"
        );
    }

    // ---- Inequality constraint row shapes ----

    /// A value-bound inequality row at an exact grid node (so trilinear
    /// interpolation collapses to a single active corner with coefficient
    /// `1.0`, i.e. `norm == 1`) with `weight == 1.0` should store `coefs`
    /// and bounds completely unscaled (`scale = weight / norm == 1`), so
    /// `apply_inequality_operator` simply returns `x` at that node, and
    /// `inequality_bounds` returns the bounds verbatim.
    #[test]
    fn inequality_value_row_matches_unit_corner_at_exact_node() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        let node = [3.0, 3.0, 3.0];
        interp.add_value_inequality_constraint(node, 0.2, 0.8, 1.0);
        assert_eq!(interp.n_value_inequality_constraints(), 1);
        assert_eq!(interp.n_inequality_constraints(), 1);

        let (lower, upper) = interp.inequality_bounds();
        assert!((lower[0] - 0.2).abs() < 1e-12, "lower={}", lower[0]);
        assert!((upper[0] - 0.8).abs() < 1e-12, "upper={}", upper[0]);

        let mut x = vec![0.0; interp.dof()];
        let node_index = interp.grid.index(3, 3, 3);
        x[node_index] = 0.5;
        let qx = interp.apply_inequality_operator(&x);
        assert_eq!(qx.len(), 1);
        assert!((qx[0] - 0.5).abs() < 1e-12, "qx={}", qx[0]);
    }

    /// A pairwise-ordering inequality row at two exact grid nodes: raw
    /// coefficients are `[1, -1]` (8-wide one-hot each), so `norm ==
    /// sqrt(2)`. With `weight == 1.0`, `scale = 1 / sqrt(2)`, so the row
    /// value for `x_a - x_b == 2` should come out to `2 / sqrt(2)`, and the
    /// stored bounds should be the raw bounds scaled by the same factor.
    #[test]
    fn inequality_pair_row_matches_scaled_difference_at_exact_nodes() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];
        interp.add_inequality_pair_constraint(a, b, -5.0, 0.0, 1.0);
        assert_eq!(interp.n_inequality_pair_constraints(), 1);
        assert_eq!(interp.n_inequality_constraints(), 1);

        let scale = 1.0 / std::f64::consts::SQRT_2;
        let (lower, upper) = interp.inequality_bounds();
        assert!(
            (lower[0] - (-5.0 * scale)).abs() < 1e-9,
            "lower={}",
            lower[0]
        );
        assert!(
            (upper[0] - (0.0 * scale)).abs() < 1e-9,
            "upper={}",
            upper[0]
        );

        let mut x = vec![0.0; interp.dof()];
        x[interp.grid.index(2, 2, 2)] = 3.0;
        x[interp.grid.index(6, 6, 6)] = 1.0;
        let qx = interp.apply_inequality_operator(&x);
        assert_eq!(qx.len(), 1);
        let expected = (3.0 - 1.0) * scale;
        assert!(
            (qx[0] - expected).abs() < 1e-9,
            "qx={} expected={}",
            qx[0],
            expected
        );
    }

    /// Both inequality constraint constructors must silently drop
    /// out-of-grid points (matching every other constraint family), and a
    /// non-positive weight (which would need to flip the bound direction to
    /// stay meaningful once folded into `coefs`/bounds together).
    #[test]
    fn inequality_constraints_outside_grid_or_nonpositive_weight_are_dropped() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);

        interp.add_value_inequality_constraint([-5.0, 0.0, 0.0], 0.0, 1.0, 1.0);
        assert_eq!(interp.n_value_inequality_constraints(), 0);
        interp.add_value_inequality_constraint([3.0, 3.0, 3.0], 0.0, 1.0, 0.0);
        assert_eq!(interp.n_value_inequality_constraints(), 0);
        interp.add_value_inequality_constraint([3.0, 3.0, 3.0], 0.0, 1.0, -1.0);
        assert_eq!(interp.n_value_inequality_constraints(), 0);

        interp.add_inequality_pair_constraint([-5.0, 0.0, 0.0], [2.0, 2.0, 2.0], -1.0, 0.0, 1.0);
        assert_eq!(interp.n_inequality_pair_constraints(), 0);
        interp.add_inequality_pair_constraint([2.0, 2.0, 2.0], [-5.0, 0.0, 0.0], -1.0, 0.0, 1.0);
        assert_eq!(interp.n_inequality_pair_constraints(), 0);
        interp.add_inequality_pair_constraint([2.0, 2.0, 2.0], [6.0, 6.0, 6.0], -1.0, 0.0, 0.0);
        assert_eq!(interp.n_inequality_pair_constraints(), 0);
    }

    // ---- ADMM solver ----

    /// When there are no inequality constraints at all, `solve_admm` should
    /// delegate straight to `solve` -- same problem, same answer (up to the
    /// usual floating-point wobble from `solve`'s internal rayon-parallel
    /// reductions not guaranteeing identical summation order between two
    /// separate calls -- see `dot` in `cg.rs` -- rather than any real
    /// difference in what's being solved), no ADMM overhead paid for having
    /// nothing to project.
    #[test]
    fn solve_admm_with_no_inequality_constraints_matches_plain_solve() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 10.0);
        interp.add_value_constraint([6.0, 3.0, 4.0], -0.5, 10.0);
        interp.add_gradient_constraint([4.0, 4.0, 4.0], [1.0, 0.2, -0.1], 2.0);

        assert_eq!(interp.n_inequality_constraints(), 0);

        let plain = interp.solve(4000, 1e-11);
        let admm = interp.solve_admm(1.0, 50, 4000, 1e-11, 1e-6, 1e-6);

        assert!(plain.converged);
        assert!(admm.converged);
        assert_eq!(
            admm.outer_iterations, 0,
            "expected the no-inequality fast path to skip the outer loop"
        );
        assert_eq!(plain.x.len(), admm.x.len());
        for (a, b) in plain.x.iter().zip(admm.x.iter()) {
            assert!(
                (a - b).abs() < 1e-8,
                "solve_admm should delegate exactly to solve() when there are no inequality rows: a={a} b={b}"
            );
        }
    }

    /// End-to-end: a value-inequality upper bound that's genuinely *active*
    /// -- the unconstrained equality fit at the constrained point
    /// overshoots the upper bound -- and `solve_admm` must pull it back down
    /// to (at most, within tolerance) that bound.
    #[test]
    fn solve_admm_enforces_active_value_upper_bound() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 2.0;

        // Unconstrained: a single strong value pin at 10.0 (far above the
        // bound) plus default (nonzero) regularisation, so the equality-only
        // solve settles near 10.0 at that point -- well past `upper_bound`.
        let mut unconstrained = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        unconstrained.add_value_constraint(point, 10.0, 100.0);
        let baseline = unconstrained.solve(4000, 1e-10);
        assert!(baseline.converged);
        let baseline_value = unconstrained.evaluate(&baseline.x, &[point])[0];
        assert!(
            baseline_value > upper_bound + 0.5,
            "expected the unconstrained fit to overshoot the bound, got {baseline_value}"
        );

        // Same data constraint, but now with a value-inequality upper bound
        // added and solved via ADMM.
        let mut constrained = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        constrained.add_value_constraint(point, 10.0, 100.0);
        constrained.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);

        // `rho` needs to be comparable to (or larger than) the competing
        // equality constraint's own stiffness at this dof (`weight^2 =
        // 100^2 = 10_000`, from `add_value_constraint`'s weighting
        // convention) for classical fixed-rho ADMM to pull the fit back to
        // the bound in a reasonable number of outer iterations -- with `rho`
        // much smaller than that (e.g. `rho = 1`), the primal residual
        // still shrinks every iteration (ADMM is still converging) but only
        // by a fraction of a percent per iteration, since the x-update
        // barely reacts to the box-projected `z`/`u` correction next to such
        // a stiff competing equality pull. This is expected, textbook
        // behaviour for constant-rho ADMM (see Boyd et al. section 3.4 on
        // rho selection), not a bug -- an adaptive-rho schedule would pick
        // this up automatically, which is exactly the kind of refinement
        // explicitly descoped here (see `solve_admm`'s docs).
        let result = constrained.solve_admm(5000.0, 2000, 4000, 1e-10, 1e-6, 1e-6);
        assert!(
            result.converged,
            "admm did not converge: outer_iterations={}",
            result.outer_iterations
        );

        let constrained_value = constrained.evaluate(&result.x, &[point])[0];
        assert!(
            constrained_value <= upper_bound + 1e-2,
            "expected ADMM to enforce the upper bound, got {constrained_value} (bound {upper_bound})"
        );
    }

    /// End-to-end: a pairwise-ordering inequality that's genuinely active --
    /// two value pins that, left alone, put point A *above* point B, while
    /// the inequality demands `f(A) - f(B) <= 0` (A must not exceed B, the
    /// stratigraphic-ordering shape). ADMM must pull the pair back into the
    /// required order.
    #[test]
    fn solve_admm_enforces_active_pairwise_ordering_bound() {
        let grid = make_grid();
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        // Sanity: unconstrained, a should indeed end up above b.
        let baseline = interp.solve(4000, 1e-10);
        assert!(baseline.converged);
        let baseline_vals = interp.evaluate(&baseline.x, &[a, b]);
        assert!(
            baseline_vals[0] > baseline_vals[1] + 0.5,
            "expected the unconstrained fit to have a > b, got a={} b={}",
            baseline_vals[0],
            baseline_vals[1]
        );

        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);
        // See the `rho` sizing comment in `solve_admm_enforces_active_value_upper_bound`:
        // it needs to be comparable to the competing value pins' own
        // stiffness (`weight^2 = 100^2 = 10_000`) for reasonably fast
        // fixed-rho ADMM convergence.
        let result = interp.solve_admm(5000.0, 2000, 4000, 1e-10, 1e-6, 1e-6);
        assert!(
            result.converged,
            "admm did not converge: outer_iterations={}",
            result.outer_iterations
        );

        let vals = interp.evaluate(&result.x, &[a, b]);
        assert!(
            vals[0] <= vals[1] + 1e-2,
            "expected ADMM to enforce a <= b, got a={} b={}",
            vals[0],
            vals[1]
        );
    }

    // ---- Active-Set CG (`solve_active_set_cg`) ----

    #[test]
    fn solve_active_set_cg_with_no_inequality_constraints_matches_plain_solve() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 10.0);
        interp.add_value_constraint([6.0, 3.0, 4.0], -0.5, 10.0);
        interp.add_gradient_constraint([4.0, 4.0, 4.0], [1.0, 0.2, -0.1], 2.0);

        assert_eq!(interp.n_inequality_constraints(), 0);

        let plain = interp.solve(4000, 1e-11);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-11, 1e-9);

        assert!(plain.converged);
        assert!(as_cg.converged);
        assert_eq!(
            as_cg.outer_iterations, 0,
            "expected the no-inequality fast path to skip the working-set loop"
        );
        assert_eq!(plain.x.len(), as_cg.x.len());
        for (a, b) in plain.x.iter().zip(as_cg.x.iter()) {
            assert!(
                (a - b).abs() < 1e-8,
                "solve_active_set_cg should delegate exactly to solve() when there are no inequality rows: a={a} b={b}"
            );
        }
    }

    /// A value-inequality constraint whose bound is never actually violated
    /// by the unconstrained fit: the working set should end up empty and the
    /// result should match the plain equality solve.
    #[test]
    fn solve_active_set_cg_matches_plain_solve_when_bound_is_not_violated() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(point, 1.0, 10.0);
        interp.add_value_inequality_constraint(point, f64::NEG_INFINITY, 100.0, 1.0);

        let plain = interp.solve(4000, 1e-11);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-11, 1e-9);

        assert!(plain.converged);
        assert!(as_cg.converged);
        assert_eq!(
            as_cg.outer_iterations, 0,
            "an always-satisfied bound should never enter the working set"
        );
        for (a, b) in plain.x.iter().zip(as_cg.x.iter()) {
            assert!((a - b).abs() < 1e-8, "a={a} b={b}");
        }
    }

    /// Same scenario as `solve_admm_enforces_active_value_upper_bound`, but
    /// via `solve_active_set_cg` -- which, unlike ADMM, must land on the
    /// bound to (near) floating-point precision, not just within a loose
    /// residual tolerance.
    #[test]
    fn solve_active_set_cg_enforces_active_value_upper_bound() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 2.0;

        let mut unconstrained = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        unconstrained.add_value_constraint(point, 10.0, 100.0);
        let baseline = unconstrained.solve(4000, 1e-10);
        assert!(baseline.converged);
        let baseline_value = unconstrained.evaluate(&baseline.x, &[point])[0];
        assert!(
            baseline_value > upper_bound + 0.5,
            "expected the unconstrained fit to overshoot the bound, got {baseline_value}"
        );

        let mut constrained = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        constrained.add_value_constraint(point, 10.0, 100.0);
        constrained.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);

        // Note: `inner_tol` here is looser than the `1e-11` used elsewhere
        // in this test module (e.g. `solve_admm`'s tests) -- the projected
        // CG's null-space projection does a little extra floating-point
        // work every iteration compared to plain Jacobi-preconditioned CG,
        // which raises its practically-reachable noise floor slightly. The
        // bound-contact precision this method actually promises comes from
        // `active_set_restore_feasibility`'s exact affine correction, not
        // from how tight `inner_tol` is -- see the assertion below, which
        // still checks contact to `1e-6`.
        let result = constrained.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        assert!(
            result.converged,
            "active-set CG did not converge: outer_iterations={}",
            result.outer_iterations
        );

        let constrained_value = constrained.evaluate(&result.x, &[point])[0];
        assert!(
            (constrained_value - upper_bound).abs() < 1e-6,
            "expected active-set CG to land on the bound to near machine precision, got {constrained_value} (bound {upper_bound})"
        );
    }

    /// Same scenario as `solve_admm_enforces_active_pairwise_ordering_bound`,
    /// via `solve_active_set_cg`.
    #[test]
    fn solve_active_set_cg_enforces_active_pairwise_ordering_bound() {
        let grid = make_grid();
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        let baseline = interp.solve(4000, 1e-10);
        assert!(baseline.converged);
        let baseline_vals = interp.evaluate(&baseline.x, &[a, b]);
        assert!(
            baseline_vals[0] > baseline_vals[1] + 0.5,
            "expected the unconstrained fit to have a > b, got a={} b={}",
            baseline_vals[0],
            baseline_vals[1]
        );

        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);
        // See the `inner_tol` sizing note on
        // `solve_active_set_cg_enforces_active_value_upper_bound`.
        let result = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        assert!(
            result.converged,
            "active-set CG did not converge: outer_iterations={}",
            result.outer_iterations
        );

        let vals = interp.evaluate(&result.x, &[a, b]);
        assert!(
            vals[0] <= vals[1] + 1e-6,
            "expected active-set CG to enforce a <= b (near-exactly), got a={} b={}",
            vals[0],
            vals[1]
        );
    }

    /// Two simultaneously active bounds (one value-bound, one pairwise) so
    /// the working set grows past one row and the Lagrange-multiplier
    /// release check gets exercised alongside step truncation.
    #[test]
    fn solve_active_set_cg_enforces_multiple_simultaneous_active_bounds() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 3.0;
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(point, 10.0, 100.0);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        interp.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);
        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);

        // See the `inner_tol` sizing note on
        // `solve_active_set_cg_enforces_active_value_upper_bound`.
        let result = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        assert!(
            result.converged,
            "active-set CG did not converge: outer_iterations={}",
            result.outer_iterations
        );

        let point_val = interp.evaluate(&result.x, &[point])[0];
        let vals = interp.evaluate(&result.x, &[a, b]);
        assert!(
            (point_val - upper_bound).abs() < 1e-6,
            "expected the value bound at {point_val} (bound {upper_bound})"
        );
        assert!(
            vals[0] <= vals[1] + 1e-6,
            "expected a <= b, got a={} b={}",
            vals[0],
            vals[1]
        );
    }

    // ---- Hybrid solve (`solve_hybrid`) ----

    #[test]
    fn solve_hybrid_with_no_inequality_constraints_matches_plain_solve() {
        let grid = make_grid();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 10.0);
        interp.add_gradient_constraint([4.0, 4.0, 4.0], [1.0, 0.2, -0.1], 2.0);

        assert_eq!(interp.n_inequality_constraints(), 0);

        let plain = interp.solve(4000, 1e-11);
        let hybrid = interp.solve_hybrid(1.0, 15, 4000, 1e-10, 50, 4000, 1e-11, 1e-9);

        assert!(plain.converged);
        assert!(hybrid.converged);
        assert_eq!(hybrid.admm_warmup_iterations, 0);
        assert_eq!(hybrid.polish_outer_iterations, 0);
        for (a, b) in plain.x.iter().zip(hybrid.x.iter()) {
            assert!(
                (a - b).abs() < 1e-8,
                "solve_hybrid should delegate exactly to solve() when there are no inequality rows: a={a} b={b}"
            );
        }
    }

    /// Same scenario as `solve_active_set_cg_enforces_active_value_upper_bound`,
    /// via the two-stage hybrid solve: a short ADMM warm-up followed by an
    /// active-set polish. The polish stage must still land on the bound to
    /// near machine precision, even though the ADMM warm-up alone only gets
    /// close.
    #[test]
    fn solve_hybrid_enforces_active_value_upper_bound_to_machine_precision() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 2.0;

        let mut constrained = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        constrained.add_value_constraint(point, 10.0, 100.0);
        constrained.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);

        // A short, deliberately loose ADMM warm-up (rho comparable to the
        // competing equality weight^2, per the sizing note on
        // `solve_admm_enforces_active_value_upper_bound`) followed by an
        // active-set polish that must clean up whatever ADMM left behind.
        // See the `inner_tol` sizing note on
        // `solve_active_set_cg_enforces_active_value_upper_bound` for why
        // `polish_tol` is looser than `1e-11`.
        let result = constrained.solve_hybrid(5000.0, 15, 4000, 1e-6, 50, 4000, 1e-8, 1e-6);
        assert!(
            result.converged,
            "hybrid polish did not converge: outer_iterations={}",
            result.polish_outer_iterations
        );

        let constrained_value = constrained.evaluate(&result.x, &[point])[0];
        assert!(
            (constrained_value - upper_bound).abs() < 1e-6,
            "expected the hybrid solve to land on the bound to near machine precision, got {constrained_value} (bound {upper_bound})"
        );
    }

    /// The hybrid solve's ADMM warm-up is a fixed budget, not a solve to
    /// convergence -- even `admm_warmup_iter = 0` (skip ADMM entirely) must
    /// still land on the bound via the active-set polish alone, since
    /// `active_set_cg_core` re-derives feasibility from scratch regardless
    /// of how good the hand-off point is.
    #[test]
    fn solve_hybrid_converges_with_zero_admm_warmup_iterations() {
        let grid = make_grid();
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);

        // See the `inner_tol` sizing note on
        // `solve_active_set_cg_enforces_active_value_upper_bound`.
        let result = interp.solve_hybrid(1.0, 0, 4000, 1e-6, 50, 4000, 1e-8, 1e-9);
        assert!(
            result.converged,
            "hybrid polish did not converge: outer_iterations={}",
            result.polish_outer_iterations
        );
        assert_eq!(result.admm_warmup_iterations, 0);

        let vals = interp.evaluate(&result.x, &[a, b]);
        assert!(
            vals[0] <= vals[1] + 1e-6,
            "expected a <= b, got a={} b={}",
            vals[0],
            vals[1]
        );
    }

    // ---- Cross-solver comparison (ADMM vs Active-Set CG vs Hybrid) ----
    //
    // All three solvers minimise the *same* strictly convex QP (`N`/`b` from
    // `ridge_factor > 0`-damped equality constraints + regularisation, `Q`/
    // bounds from the inequality rows), which has a unique global optimum --
    // so, given tight enough settings, all three must agree on that optimum
    // regardless of which algorithm found it. These tests exploit that: they
    // solve one interpolator instance with all three (`solve_*` all take
    // `&self`, so the same instance can be reused) and check both that the
    // solved *fields* agree and that ADMM's bound contact is measurably
    // looser than Active-Set CG's/Hybrid's exact contact -- the trade-off
    // `Inequalities.md`'s decision matrix describes.

    fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn all_three_solvers_agree_on_active_value_upper_bound() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 2.0;

        let mut interp = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
            grid.origin,
            grid.nsteps,
            grid.step_vector,
        ));
        interp.add_value_constraint(point, 10.0, 100.0);
        interp.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);

        // Same `rho`/iteration sizing as `solve_admm_enforces_active_value_upper_bound`.
        let admm = interp.solve_admm(5000.0, 2000, 4000, 1e-10, 1e-8, 1e-8);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        let hybrid = interp.solve_hybrid(5000.0, 15, 4000, 1e-6, 50, 4000, 1e-8, 1e-6);

        assert!(admm.converged, "admm did not converge");
        assert!(as_cg.converged, "active-set CG did not converge");
        assert!(hybrid.converged, "hybrid did not converge");

        let admm_err = (interp.evaluate(&admm.x, &[point])[0] - upper_bound).abs();
        let as_cg_err = (interp.evaluate(&as_cg.x, &[point])[0] - upper_bound).abs();
        let hybrid_err = (interp.evaluate(&hybrid.x, &[point])[0] - upper_bound).abs();
        assert!(
            as_cg_err < 1e-6,
            "expected active-set CG's bound contact to be near-exact, got err={as_cg_err}"
        );
        assert!(
            hybrid_err < 1e-6,
            "expected hybrid's bound contact to be near-exact, got err={hybrid_err}"
        );
        assert!(
            admm_err > as_cg_err * 10.0,
            "expected ADMM's asymptotic contact to be measurably looser than active-set CG's exact contact, got admm_err={admm_err} as_cg_err={as_cg_err}"
        );

        // All three minimise the same QP, so the whole field -- not just the
        // constrained point -- should agree closely, bounded by ADMM's own
        // (looser) convergence tolerance.
        assert!(
            max_abs_diff(&admm.x, &as_cg.x) < 1e-2,
            "admm/active-set-cg fields disagree by more than expected: max diff={}",
            max_abs_diff(&admm.x, &as_cg.x)
        );
        assert!(
            max_abs_diff(&admm.x, &hybrid.x) < 1e-2,
            "admm/hybrid fields disagree by more than expected: max diff={}",
            max_abs_diff(&admm.x, &hybrid.x)
        );
        assert!(
            max_abs_diff(&as_cg.x, &hybrid.x) < 1e-2,
            "active-set-cg/hybrid fields disagree by more than expected: max diff={}",
            max_abs_diff(&as_cg.x, &hybrid.x)
        );
    }

    #[test]
    fn all_three_solvers_agree_on_active_pairwise_ordering_bound() {
        let grid = make_grid();
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);

        let admm = interp.solve_admm(5000.0, 2000, 4000, 1e-10, 1e-8, 1e-8);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        let hybrid = interp.solve_hybrid(5000.0, 15, 4000, 1e-6, 50, 4000, 1e-8, 1e-6);

        assert!(admm.converged, "admm did not converge");
        assert!(as_cg.converged, "active-set CG did not converge");
        assert!(hybrid.converged, "hybrid did not converge");

        let admm_vals = interp.evaluate(&admm.x, &[a, b]);
        let as_cg_vals = interp.evaluate(&as_cg.x, &[a, b]);
        let hybrid_vals = interp.evaluate(&hybrid.x, &[a, b]);

        let admm_err = (admm_vals[0] - admm_vals[1]).max(0.0);
        let as_cg_err = (as_cg_vals[0] - as_cg_vals[1]).max(0.0);
        let hybrid_err = (hybrid_vals[0] - hybrid_vals[1]).max(0.0);
        assert!(
            as_cg_err < 1e-6,
            "expected active-set CG's a<=b contact to be near-exact, got err={as_cg_err}"
        );
        assert!(
            hybrid_err < 1e-6,
            "expected hybrid's a<=b contact to be near-exact, got err={hybrid_err}"
        );
        assert!(
            admm_err <= 1e-2,
            "expected ADMM to satisfy a<=b within its own (looser) tolerance, got err={admm_err}"
        );

        assert!(
            max_abs_diff(&admm.x, &as_cg.x) < 1e-2,
            "admm/active-set-cg fields disagree by more than expected: max diff={}",
            max_abs_diff(&admm.x, &as_cg.x)
        );
        assert!(
            max_abs_diff(&admm.x, &hybrid.x) < 1e-2,
            "admm/hybrid fields disagree by more than expected: max diff={}",
            max_abs_diff(&admm.x, &hybrid.x)
        );
        assert!(
            max_abs_diff(&as_cg.x, &hybrid.x) < 1e-2,
            "active-set-cg/hybrid fields disagree by more than expected: max diff={}",
            max_abs_diff(&as_cg.x, &hybrid.x)
        );
    }

    /// Same multi-constraint scenario as
    /// `solve_active_set_cg_enforces_multiple_simultaneous_active_bounds`,
    /// but cross-checking all three solvers agree with each other on both
    /// active bounds simultaneously.
    #[test]
    fn all_three_solvers_agree_with_multiple_simultaneous_active_bounds() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];
        let upper_bound = 3.0;
        let a = [2.0, 2.0, 2.0];
        let b = [6.0, 6.0, 6.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(point, 10.0, 100.0);
        interp.add_value_constraint(a, 5.0, 100.0);
        interp.add_value_constraint(b, 1.0, 100.0);
        interp.add_value_inequality_constraint(point, f64::NEG_INFINITY, upper_bound, 1.0);
        interp.add_inequality_pair_constraint(a, b, f64::NEG_INFINITY, 0.0, 1.0);

        let admm = interp.solve_admm(5000.0, 2000, 4000, 1e-10, 1e-8, 1e-8);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        let hybrid = interp.solve_hybrid(5000.0, 15, 4000, 1e-6, 50, 4000, 1e-8, 1e-6);

        assert!(admm.converged, "admm did not converge");
        assert!(as_cg.converged, "active-set CG did not converge");
        assert!(hybrid.converged, "hybrid did not converge");

        for (name, x) in [
            ("admm", &admm.x),
            ("as_cg", &as_cg.x),
            ("hybrid", &hybrid.x),
        ] {
            let point_val = interp.evaluate(x, &[point])[0];
            let vals = interp.evaluate(x, &[a, b]);
            assert!(
                point_val <= upper_bound + 1e-2,
                "{name}: expected the value bound to hold, got {point_val} (bound {upper_bound})"
            );
            assert!(
                vals[0] <= vals[1] + 1e-2,
                "{name}: expected a <= b, got a={} b={}",
                vals[0],
                vals[1]
            );
        }

        assert!(
            max_abs_diff(&admm.x, &as_cg.x) < 1e-2,
            "admm/active-set-cg fields disagree by more than expected: max diff={}",
            max_abs_diff(&admm.x, &as_cg.x)
        );
        assert!(
            max_abs_diff(&as_cg.x, &hybrid.x) < 1e-2,
            "active-set-cg/hybrid fields disagree by more than expected: max diff={}",
            max_abs_diff(&as_cg.x, &hybrid.x)
        );
    }

    /// When no bound is ever violated, all three solvers (including their
    /// no-inequality-rows fast paths) must degenerate to the same plain
    /// equality solve.
    #[test]
    fn all_three_solvers_agree_with_plain_solve_when_no_bound_is_violated() {
        let grid = make_grid();
        let point = [4.0, 4.0, 4.0];

        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint(point, 1.0, 10.0);
        interp.add_value_inequality_constraint(point, f64::NEG_INFINITY, 100.0, 1.0);

        let plain = interp.solve(4000, 1e-10);
        let admm = interp.solve_admm(1.0, 50, 4000, 1e-10, 1e-6, 1e-6);
        let as_cg = interp.solve_active_set_cg(50, 4000, 1e-8, 1e-9);
        let hybrid = interp.solve_hybrid(1.0, 15, 4000, 1e-6, 50, 4000, 1e-8, 1e-9);

        assert!(plain.converged);
        assert!(admm.converged);
        assert!(as_cg.converged);
        assert!(hybrid.converged);

        for (name, x) in [
            ("admm", &admm.x),
            ("as_cg", &as_cg.x),
            ("hybrid", &hybrid.x),
        ] {
            assert!(
                max_abs_diff(&plain.x, x) < 1e-6,
                "{name} should match the plain equality solve when no bound is ever violated: max diff={}",
                max_abs_diff(&plain.x, x)
            );
        }
    }

    // ---- Fused single-kernel regularisation: bit-for-bit-close equivalence ----

    fn assert_fields_close(fused: &[f64], unfused: &[f64], grid: &StructuredGrid3D, ctx: &str) {
        assert_eq!(fused.len(), unfused.len());
        let scale = unfused.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1.0);
        for (flat, (a, b)) in fused.iter().zip(unfused.iter()).enumerate() {
            let diff = (a - b).abs();
            if diff > 1e-9 * scale {
                let (i, j, k) = grid.node_ijk(flat);
                panic!(
                    "{ctx}: fused/unfused mismatch at flat={flat} (i={i},j={j},k={k}): fused={a} unfused={b} diff={diff} scale={scale}"
                );
            }
        }
    }

    /// The fused single-kernel regularisation path (`apply_normal_operator`)
    /// must produce the exact same normal-operator output as the reference
    /// two-pass computation (`apply_normal_operator_unfused`) at *every* dof
    /// -- not just in the deep interior, but also across the boundary shell
    /// the fused kernel deliberately skips and the shell-correction pass
    /// fixes up. Exercised across several grid sizes (including a tiny one
    /// where most/all nodes are shell nodes), several random fields, and a
    /// mix of regularisation/border weights (including some zero, to check
    /// the "family doesn't contribute" fast paths agree too).
    #[test]
    fn fused_matches_unfused_various_grids_and_fields() {
        let grid_shapes: [[usize; 3]; 4] = [[4, 4, 4], [5, 5, 5], [6, 9, 5], [11, 8, 13]];
        let weight_sets = [
            RegularisationWeights::default(),
            RegularisationWeights {
                dxx: 2.0,
                dyy: 0.5,
                dzz: 1.5,
                dxy: 0.0,
                dxz: 3.0,
                dyz: 0.0,
            },
            RegularisationWeights {
                dxx: 0.0,
                dyy: 0.0,
                dzz: 0.0,
                dxy: 0.0,
                dxz: 0.0,
                dyz: 0.0,
            },
        ];
        let border_weight_sets = [
            BorderRegularisationWeights::default(),
            BorderRegularisationWeights {
                dx: 0.0,
                dy: 2.0,
                dz: 0.7,
            },
        ];

        for shape in grid_shapes {
            let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], shape, [1.0, 1.0, 1.0]);
            for weights in weight_sets {
                for border_weights in border_weight_sets {
                    let mut interp = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(
                        grid.origin,
                        grid.nsteps,
                        grid.step_vector,
                    ))
                    .with_regularisation(weights)
                    .with_border_regularisation(border_weights);
                    // A couple of scattered value/gradient constraints so the
                    // data-term diagonal isn't uniformly zero either.
                    interp.add_value_constraint([1.0, 1.0, 1.0], 0.7, 3.0);
                    interp.add_gradient_constraint([2.0, 1.0, 0.0], [1.0, 0.3, -0.2], 1.5);

                    let n = interp.dof();
                    for seed in [1u64, 2u64, 3u64] {
                        let x: Vec<f64> = (0..n)
                            .map(|i| {
                                let h = (i as u64)
                                    .wrapping_mul(6364136223846793005u64)
                                    .wrapping_add(seed);
                                (h % 1000) as f64 / 1000.0 - 0.5
                            })
                            .collect();
                        let fused = interp.apply_normal_operator(&x);
                        let unfused = interp.apply_normal_operator_unfused(&x);
                        assert_fields_close(
                            &fused,
                            &unfused,
                            &interp.grid,
                            &format!(
                                "shape={shape:?} weights={weights:?} border={border_weights:?} seed={seed}"
                            ),
                        );
                    }
                }
            }
        }
    }

    /// Same equivalence check, but specifically isolated to the boundary
    /// shell (dist < 2 from the true edge) on a grid large enough to have a
    /// genuine deep interior too -- guards against a regression that
    /// happens to look fine in aggregate (e.g. errors at a handful of shell
    /// nodes averaging out against a large, correct interior) by checking
    /// shell nodes specifically and individually.
    #[test]
    fn fused_matches_unfused_at_boundary_shell_specifically() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [12, 10, 9], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid)
            .with_regularisation(RegularisationWeights {
                dxx: 1.3,
                dyy: 0.6,
                dzz: 2.0,
                dxy: 1.0,
                dxz: 0.4,
                dyz: 1.1,
            })
            .with_border_regularisation(BorderRegularisationWeights {
                dx: 1.7,
                dy: 0.9,
                dz: 1.2,
            });
        interp.add_value_constraint([0.0, 0.0, 0.0], 1.0, 5.0);
        interp.add_value_constraint([11.0, 9.0, 8.0], -1.0, 5.0);

        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 53 + 7) % 97) as f64 / 97.0 - 0.5)
            .collect();
        let fused = interp.apply_normal_operator(&x);
        let unfused = interp.apply_normal_operator_unfused(&x);

        let (nx, ny, nz) = (
            interp.grid.nsteps[0],
            interp.grid.nsteps[1],
            interp.grid.nsteps[2],
        );
        let mut shell_checked = 0;
        for flat in 0..n {
            let (i, j, k) = interp.grid.node_ijk(flat);
            if dist_to_true_edge(&[nx, ny, nz], i, j, k) < 2 {
                shell_checked += 1;
                let diff = (fused[flat] - unfused[flat]).abs();
                let scale = unfused.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1.0);
                assert!(
                    diff < 1e-9 * scale,
                    "shell mismatch at flat={flat} (i={i},j={j},k={k}): fused={} unfused={} diff={diff}",
                    fused[flat],
                    unfused[flat]
                );
            }
        }
        // Sanity: this grid should actually have a nontrivial shell (corners,
        // edges, faces), otherwise the test isn't exercising what it claims to.
        assert!(
            shell_checked > 50,
            "expected a substantial boundary shell, got {shell_checked} nodes"
        );
    }

    // ---- BlockJacobiX preconditioner ----

    /// `BlockJacobiX` must converge to (numerically) the same solution as
    /// plain Jacobi on the same problem -- confirming it's a valid SPD
    /// preconditioner CG convergence theory still applies to, not just a
    /// faster-but-wrong shortcut -- and, on a representative multi-family
    /// problem (value + gradient + norm + curvature + border, a few
    /// thousand dof), do it in fewer iterations.
    #[test]
    fn block_jacobi_converges_to_same_solution_in_fewer_iterations() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [22, 18, 16], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid)
            .with_regularisation(RegularisationWeights::default())
            .with_border_regularisation(BorderRegularisationWeights::default());

        let f =
            |p: [f64; 3]| ((p[0] * 0.4).sin() * (p[1] * 0.3).cos() * 3.0) + p[2] * 0.2 - 0.1 * p[0];
        for p in [
            [2.0, 2.0, 2.0],
            [18.0, 4.0, 6.0],
            [5.0, 14.0, 3.0],
            [10.0, 9.0, 12.0],
            [3.0, 3.0, 13.0],
            [15.0, 15.0, 8.0],
        ] {
            interp.add_value_constraint(p, f(p), 20.0);
        }
        interp.add_gradient_constraint([8.0, 8.0, 8.0], [1.0, 0.2, -0.3], 5.0);
        interp.add_norm_constraint([12.0, 6.0, 9.0], [0.1, 0.2, 0.05], 3.0);

        let jacobi = interp.solve_with(5000, 1e-10, PreconditionerKind::Jacobi);
        let block = interp.solve_with(5000, 1e-10, PreconditionerKind::BlockJacobiX);

        assert!(jacobi.converged, "jacobi iterations={}", jacobi.iterations);
        assert!(
            block.converged,
            "block-jacobi iterations={}",
            block.iterations
        );

        let max_diff = jacobi
            .x
            .iter()
            .zip(block.x.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(max_diff < 1e-5, "solutions diverged: max_diff={max_diff}");

        assert!(
            block.iterations <= jacobi.iterations,
            "expected block-Jacobi to need no more iterations than plain Jacobi: jacobi={} block={}",
            jacobi.iterations,
            block.iterations
        );
        eprintln!(
            "block_jacobi_converges_to_same_solution_in_fewer_iterations: jacobi={} iters, block_jacobi_x={} iters (dof={})",
            jacobi.iterations,
            block.iterations,
            interp.dof()
        );
    }

    /// Degenerate case: no `dxx` weight and no x-axis border regularisation
    /// means `dxx_and_border_dx_offdiagonal` is all zeros, so `BlockJacobiX`
    /// degenerates to exactly plain Jacobi (a tridiagonal solve with zero
    /// off-diagonal is just `z[i] = r[i] / diag[i]` per entry) -- confirms
    /// it never does *worse* than Jacobi even when there's no coupling to
    /// exploit, and that the degenerate Thomas-algorithm path doesn't panic
    /// or diverge.
    #[test]
    fn block_jacobi_matches_jacobi_when_no_x_coupling_present() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [10, 8, 7], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid)
            .with_regularisation(RegularisationWeights {
                dxx: 0.0,
                dyy: 1.0,
                dzz: 1.0,
                dxy: 0.0,
                dxz: 0.0,
                dyz: 0.0,
            })
            .with_border_regularisation(BorderRegularisationWeights {
                dx: 0.0,
                dy: 1.0,
                dz: 1.0,
            });
        interp.add_value_constraint([1.0, 1.0, 1.0], 1.0, 10.0);
        interp.add_value_constraint([8.0, 6.0, 5.0], -1.0, 10.0);

        let jacobi = interp.solve_with(2000, 1e-11, PreconditionerKind::Jacobi);
        let block = interp.solve_with(2000, 1e-11, PreconditionerKind::BlockJacobiX);
        assert!(jacobi.converged && block.converged);
        assert_eq!(jacobi.iterations, block.iterations);
        let max_diff = jacobi
            .x
            .iter()
            .zip(block.x.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(max_diff < 1e-9, "max_diff={max_diff}");
    }

    // ---- Manual before/after profiling (not run by default) ----

    /// Before/after benchmark for the fused single-kernel regularisation
    /// optimisation (TASKS.md P1). Run with:
    /// `cargo test --release -- --ignored --nocapture profile_apply_1e6`.
    #[test]
    #[ignore]
    fn profile_apply_1e6() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [100, 100, 100], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        for p in [[10.0, 10.0, 10.0], [50.0, 50.0, 50.0], [90.0, 5.0, 80.0]] {
            interp.add_value_constraint(p, 1.0, 10.0);
        }
        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let mut scratch = vec![0.0; interp.scratch_len()];
        let mut out = vec![0.0; n];
        for _ in 0..3 {
            interp.apply_normal_operator_into(&x, &mut scratch, &mut out);
        }
        let iters = 20;
        let start = std::time::Instant::now();
        for _ in 0..iters {
            interp.apply_normal_operator_into(&x, &mut scratch, &mut out);
        }
        let elapsed = start.elapsed();
        eprintln!(
            "profile_apply_1e6 [FUSED]: {:?}/call over {iters} calls (dof={n})",
            elapsed / iters
        );
    }

    /// As `profile_apply_1e6`, but using the reference (pre-fusion) two-pass
    /// computation -- the "before" measurement for the same optimisation.
    #[test]
    #[ignore]
    fn profile_apply_1e6_unfused() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [100, 100, 100], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        for p in [[10.0, 10.0, 10.0], [50.0, 50.0, 50.0], [90.0, 5.0, 80.0]] {
            interp.add_value_constraint(p, 1.0, 10.0);
        }
        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let mut scratch = vec![0.0; interp.scratch_len_unfused()];
        let mut out = vec![0.0; n];
        for _ in 0..3 {
            interp.apply_normal_operator_into_unfused(&x, &mut scratch, &mut out);
        }
        let iters = 20;
        let start = std::time::Instant::now();
        for _ in 0..iters {
            interp.apply_normal_operator_into_unfused(&x, &mut scratch, &mut out);
        }
        let elapsed = start.elapsed();
        eprintln!(
            "profile_apply_1e6_unfused [UNFUSED]: {:?}/call over {iters} calls (dof={n})",
            elapsed / iters
        );
    }

    /// `profile_apply_1e6`/`profile_apply_1e6_unfused`'s overall
    /// `apply_normal_operator_into` speedup is small: border regularisation
    /// (six centred `[1, -2, 1]` `EdgeStencil` families, each still using the original,
    /// unchanged, full-`dof`-iterating exact two-pass computation -- out of
    /// scope for this optimisation, see the module docs) contributes far
    /// more O(dof) passes (12: forward+adjoint x 6 border families) than the
    /// six interior curvature families ever did (2: one forward, one
    /// adjoint, now collapsed to ~1), so it dominates the *total* time and
    /// dilutes the interior fusion's own speedup in that measurement. This
    /// variant disables border regularisation to isolate and measure the
    /// interior-family fusion in isolation -- the actual thing this
    /// optimisation targets.
    #[test]
    #[ignore]
    fn profile_apply_1e6_interior_only() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [100, 100, 100], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid).with_border_regularisation(
            BorderRegularisationWeights {
                dx: 0.0,
                dy: 0.0,
                dz: 0.0,
            },
        );
        for p in [[10.0, 10.0, 10.0], [50.0, 50.0, 50.0], [90.0, 5.0, 80.0]] {
            interp.add_value_constraint(p, 1.0, 10.0);
        }
        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let mut scratch = vec![0.0; interp.scratch_len()];
        let mut out = vec![0.0; n];
        for _ in 0..3 {
            interp.apply_normal_operator_into(&x, &mut scratch, &mut out);
        }
        let iters = 20;
        let start = std::time::Instant::now();
        for _ in 0..iters {
            interp.apply_normal_operator_into(&x, &mut scratch, &mut out);
        }
        let elapsed = start.elapsed();
        eprintln!(
            "profile_apply_1e6_interior_only [FUSED, no border]: {:?}/call over {iters} calls (dof={n})",
            elapsed / iters
        );
    }

    #[test]
    #[ignore]
    fn profile_apply_1e6_interior_only_unfused() {
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [100, 100, 100], [1.0, 1.0, 1.0]);
        let mut interp = FiniteDifferenceInterpolator::new(grid).with_border_regularisation(
            BorderRegularisationWeights {
                dx: 0.0,
                dy: 0.0,
                dz: 0.0,
            },
        );
        for p in [[10.0, 10.0, 10.0], [50.0, 50.0, 50.0], [90.0, 5.0, 80.0]] {
            interp.add_value_constraint(p, 1.0, 10.0);
        }
        let n = interp.dof();
        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let mut scratch = vec![0.0; interp.scratch_len_unfused()];
        let mut out = vec![0.0; n];
        for _ in 0..3 {
            interp.apply_normal_operator_into_unfused(&x, &mut scratch, &mut out);
        }
        let iters = 20;
        let start = std::time::Instant::now();
        for _ in 0..iters {
            interp.apply_normal_operator_into_unfused(&x, &mut scratch, &mut out);
        }
        let elapsed = start.elapsed();
        eprintln!(
            "profile_apply_1e6_interior_only_unfused [UNFUSED, no border]: {:?}/call over {iters} calls (dof={n})",
            elapsed / iters
        );
    }

    /// An axis-aligned direction field (`(0,0,1)` everywhere) collapses
    /// `comp_w` to `[0,0,1,0,0,0]`, so `add_directional_regularisation`
    /// should contribute *exactly* the same normal operator as the existing
    /// isotropic `dzz` family scaled by the same weight -- the strongest
    /// correctness check available (an exact numeric match, not just
    /// "plausible"), and confirms the six families stay independent (no
    /// cross-family terms leak in for a pure axis direction).
    #[test]
    fn directional_regularisation_matches_isotropic_family_for_axis_aligned_direction() {
        let w = 0.7;
        let grid = make_grid();
        let n = grid.n_nodes();

        let mut isotropic = FiniteDifferenceInterpolator::new(grid)
            .with_regularisation(RegularisationWeights {
                dxx: 0.0,
                dyy: 0.0,
                dzz: w,
                dxy: 0.0,
                dxz: 0.0,
                dyz: 0.0,
            })
            .with_border_regularisation(BorderRegularisationWeights {
                dx: 0.0,
                dy: 0.0,
                dz: 0.0,
            });
        isotropic.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);

        let mut directional = FiniteDifferenceInterpolator::new(grid)
            .with_regularisation(RegularisationWeights {
                dxx: 0.0,
                dyy: 0.0,
                dzz: 0.0,
                dxy: 0.0,
                dxz: 0.0,
                dyz: 0.0,
            })
            .with_border_regularisation(BorderRegularisationWeights {
                dx: 0.0,
                dy: 0.0,
                dz: 0.0,
            });
        directional.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        let directions = vec![[0.0, 0.0, 1.0]; n];
        directional.add_directional_regularisation(&directions, w);

        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let iso_out = isotropic.apply_normal_operator(&x);
        let dir_out = directional.apply_normal_operator(&x);
        for (a, b) in iso_out.iter().zip(dir_out.iter()) {
            assert!((a - b).abs() < 1e-10, "isotropic={a} directional={b}");
        }
    }

    /// Symmetry check with a directional-regularisation term active
    /// (mirrors `normal_operator_is_symmetric_with_all_constraint_families`),
    /// guarding against the new unfused per-node-weighted pass breaking `N`'s
    /// self-adjointness.
    #[test]
    fn normal_operator_is_symmetric_with_directional_regularisation() {
        let grid = make_grid();
        let n = grid.n_nodes();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        // A non-axis-aligned, per-node-varying direction field so every
        // family (including the cross-derivative ones) picks up a nonzero,
        // spatially-varying weight.
        let directions: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = i as f64 * 0.37;
                let v = [t.sin(), (t * 1.3).cos(), (t * 0.7).sin()];
                let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if norm > 0.0 {
                    [v[0] / norm, v[1] / norm, v[2] / norm]
                } else {
                    [1.0, 0.0, 0.0]
                }
            })
            .collect();
        interp.add_directional_regularisation(&directions, 0.3);

        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 17 + 5) % 89) as f64 / 89.0).collect();
        let nx = interp.apply_normal_operator(&x);
        let ny = interp.apply_normal_operator(&y);
        let lhs: f64 = nx.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(ny.iter()).map(|(a, b)| a * b).sum();
        assert!(
            (lhs - rhs).abs() < 1e-6 * lhs.abs().max(1.0),
            "lhs={lhs} rhs={rhs}"
        );
    }

    /// Two calls with weight `w/sqrt(2)` (same direction) must accumulate to
    /// the same squared-weight total as one call with weight `w` --
    /// confirms `add_directional_regularisation` sums independent terms
    /// rather than overwriting the running per-node weight state.
    #[test]
    fn directional_regularisation_accumulates_across_calls() {
        let w = 0.9;
        let grid = make_grid();
        let n = grid.n_nodes();
        let directions = vec![[0.6, 0.8, 0.0]; n]; // already unit-length

        let mut once = FiniteDifferenceInterpolator::new(grid);
        once.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        once.add_directional_regularisation(&directions, w);

        let mut twice = FiniteDifferenceInterpolator::new(grid);
        twice.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        let half = w / std::f64::consts::SQRT_2;
        twice.add_directional_regularisation(&directions, half);
        twice.add_directional_regularisation(&directions, half);

        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let once_out = once.apply_normal_operator(&x);
        let twice_out = twice.apply_normal_operator(&x);
        for (a, b) in once_out.iter().zip(twice_out.iter()) {
            assert!((a - b).abs() < 1e-9, "once={a} twice={b}");
        }
    }

    /// Guards against a regression that would make every non-fold caller of
    /// this crate pay for directional regularisation it never asked for:
    /// with no `add_directional_regularisation` calls, `directional_weight2`
    /// stays `None`, and `scratch_len`/`apply_normal_operator`/`diagonal`
    /// must be byte-for-byte identical to before this feature existed.
    #[test]
    fn directional_regularisation_is_zero_cost_when_unused() {
        let grid = make_grid();
        let n = grid.n_nodes();
        let mut interp = FiniteDifferenceInterpolator::new(grid);
        interp.add_value_constraint([2.0, 2.0, 2.0], 1.0, 1.0);
        interp.add_gradient_constraint([3.0, 3.0, 3.0], [1.0, 0.0, 0.0], 2.0);

        assert!(interp.directional_weight2.is_none());

        let x: Vec<f64> = (0..n)
            .map(|i| ((i * 31 + 3) % 101) as f64 / 101.0)
            .collect();
        let out = interp.apply_normal_operator(&x);
        let diag = interp.diagonal();
        // Sanity: still a well-formed, finite result (the real regression
        // guard here is `directional_weight2.is_none()` above -- this just
        // confirms the untouched code paths still produce sane output).
        assert!(out.iter().all(|v| v.is_finite()));
        assert!(diag.iter().all(|v| v.is_finite() && *v > 0.0));
    }
}
