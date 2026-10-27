//! Minimal matrix-free (preconditioned) conjugate gradient solver.
//!
//! `apply`/`precondition` write their result into a caller-supplied output
//! buffer instead of returning a freshly allocated `Vec`, and the solver
//! itself allocates its five `dof`-length work vectors (`x`, `r`, `z`, `p`,
//! `ap`) exactly once regardless of `max_iter`. For a problem of this shape
//! (many cheap iterations over a modest-sized dense array) allocator churn
//! -- not arithmetic -- is the dominant cost if each iteration allocates
//! fresh vectors, so this is the difference between a solver that's
//! actually fast and one that merely looks matrix-free on paper.

use rayon::prelude::*;

use crate::PAR_MIN_LEN;

pub struct CgResult {
    pub x: Vec<f64>,
    pub iterations: usize,
    pub residual_norm: f64,
    pub converged: bool,
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter()
        .zip(b.par_iter())
        .with_min_len(PAR_MIN_LEN)
        .map(|(x, y)| x * y)
        .sum()
}

/// Solve `N x = rhs` with plain (unpreconditioned) conjugate gradient, `N`
/// given implicitly by `apply(x, out)` computing `out = N x`.
pub fn conjugate_gradient<F>(
    apply: F,
    rhs: &[f64],
    x0: Option<Vec<f64>>,
    max_iter: usize,
    tol: f64,
) -> CgResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    conjugate_gradient_preconditioned(
        apply,
        |r: &[f64], out: &mut [f64]| out.copy_from_slice(r),
        rhs,
        x0,
        max_iter,
        tol,
    )
}

/// Solve `N x = rhs` with preconditioned conjugate gradient: `N` is given
/// implicitly by `apply(x, out)` (`out = N x`), and the preconditioner
/// `M^-1` implicitly by `precondition(r, out)` (`out = M^-1 r`), e.g. a
/// Jacobi/diagonal solve (see `FiniteDifferenceInterpolator::diagonal`).
/// Passing `|r, out| out.copy_from_slice(r)` recovers plain CG.
pub fn conjugate_gradient_preconditioned<F, P>(
    mut apply: F,
    mut precondition: P,
    rhs: &[f64],
    x0: Option<Vec<f64>>,
    max_iter: usize,
    tol: f64,
) -> CgResult
where
    F: FnMut(&[f64], &mut [f64]),
    P: FnMut(&[f64], &mut [f64]),
{
    let n = rhs.len();
    let mut x = x0.unwrap_or_else(|| vec![0.0; n]);
    let mut r = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut ap = vec![0.0; n];

    apply(&x, &mut ap); // ap = N x0
    r.par_iter_mut()
        .zip(rhs.par_iter())
        .zip(ap.par_iter())
        .with_min_len(PAR_MIN_LEN)
        .for_each(|((ri, b), a)| *ri = b - a);
    precondition(&r, &mut z);
    p.copy_from_slice(&z);
    let mut rz_old = dot(&r, &z);

    let rhs_norm = dot(rhs, rhs).sqrt().max(1e-300);
    let mut iterations = 0usize;
    let mut converged = dot(&r, &r).sqrt() / rhs_norm <= tol;

    while !converged && iterations < max_iter {
        apply(&p, &mut ap); // ap = N p
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-300 {
            break;
        }
        let alpha = rz_old / p_ap;
        x.par_iter_mut()
            .zip(p.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(xi, pi)| *xi += alpha * pi);
        r.par_iter_mut()
            .zip(ap.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(ri, api)| *ri -= alpha * api);

        iterations += 1;
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / rhs_norm <= tol {
            converged = true;
            break;
        }
        precondition(&r, &mut z);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        p.par_iter_mut()
            .zip(z.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(pi, zi)| *pi = zi + beta * *pi);
        rz_old = rz_new;
    }

    CgResult {
        x,
        iterations,
        residual_norm: dot(&r, &r).sqrt(),
        converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_diagonal_system() {
        // N = diag(1,2,3,4), rhs = [1,2,3,4] => x = [1,1,1,1]
        let diag = [1.0, 2.0, 3.0, 4.0];
        let apply = |x: &[f64], out: &mut [f64]| {
            for (o, (xi, d)) in out.iter_mut().zip(x.iter().zip(diag.iter())) {
                *o = xi * d;
            }
        };
        let rhs = vec![1.0, 2.0, 3.0, 4.0];
        let result = conjugate_gradient(apply, &rhs, None, 100, 1e-10);
        assert!(result.converged);
        for v in result.x {
            assert!((v - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn jacobi_preconditioner_matches_plain_cg_solution() {
        let diag = [1.0, 2.0, 3.0, 4.0, 5.0];
        let apply = |x: &[f64], out: &mut [f64]| {
            for (o, (xi, d)) in out.iter_mut().zip(x.iter().zip(diag.iter())) {
                *o = xi * d;
            }
        };
        let rhs = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        let precondition = |r: &[f64], out: &mut [f64]| {
            for (o, (ri, d)) in out.iter_mut().zip(r.iter().zip(diag.iter())) {
                *o = ri / d;
            }
        };

        let plain = conjugate_gradient(apply, &rhs, None, 100, 1e-12);
        let preconditioned =
            conjugate_gradient_preconditioned(apply, precondition, &rhs, None, 100, 1e-12);

        assert!(preconditioned.iterations <= plain.iterations);
        for (a, b) in plain.x.iter().zip(preconditioned.x.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }
}
