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

use crate::error::CoreError;
use crate::memory::MemTracker;
use crate::par::*;
use crate::scalar::Scalar;
use std::any::Any;

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
    let mut p = vec![0.0; n];
    // `ap` also holds `z`: `ap` is free between the update of `r` and the
    // next operator call, so no fifth work vector is needed.
    let mut ap = vec![0.0; n];

    apply(&x, &mut ap); // ap = N x0
    r.par_iter_mut()
        .zip(rhs.par_iter())
        .zip(ap.par_iter())
        .with_min_len(PAR_MIN_LEN)
        .for_each(|((ri, b), a)| *ri = b - a);
    precondition(&r, &mut ap); // ap is z here
    p.copy_from_slice(&ap);
    let mut rz_old = dot(&r, &ap);

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
        precondition(&r, &mut ap); // ap is z here
        let rz_new = dot(&r, &ap);
        let beta = rz_new / rz_old;
        p.par_iter_mut()
            .zip(ap.par_iter())
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

/// State of a CG run that is done: why it stopped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CgStatus {
    /// More iterations are possible. `residual` is relative to `|rhs|`.
    Running { iteration: usize, residual: f64 },
    /// The relative residual is at or below the tolerance.
    Converged,
    /// The total iteration limit is reached.
    MaxIter,
    /// `p.Ap` is zero, or (single precision only) the residual is NaN or
    /// does not fall any more. The solve cannot go on.
    Breakdown,
}

/// Below this relative residual the single precision stagnation test runs.
const F32_STAGNATION_LEVEL: f64 = 1e-4;
/// Number of iterations without a 1 percent fall that means stagnation.
const F32_STAGNATION_ITERS: usize = 50;

/// `true` if the vector element type `S` is `f32`.
#[inline]
pub(crate) fn is_single<S: Scalar>() -> bool {
    std::mem::size_of::<S>() == 4
}

/// `sum(a * b)`, added in `f64`.
#[inline]
pub(crate) fn dot_s<S: Scalar>(a: &[S], b: &[S]) -> f64 {
    a.par_iter()
        .zip(b.par_iter())
        .with_min_len(PAR_MIN_LEN)
        .map(|(x, y)| x.to_f64() * y.to_f64())
        .sum::<f64>()
}

/// Make the Jacobi factors `1 / diag` of type `S`. A diagonal value of zero
/// or less gives the factor 1 (no scaling). For `S = f64` the input vector is
/// changed in place and returned (no copy). For `S = f32` the function makes
/// a new `f32` vector from `tracker` and gives the `f64` bytes back.
pub(crate) fn jacobi_inverse<S: Scalar>(
    tracker: &mut MemTracker,
    mut diag: Vec<f64>,
) -> Result<Vec<S>, CoreError> {
    diag.par_iter_mut()
        .with_min_len(PAR_MIN_LEN)
        .for_each(|d| *d = if *d > 0.0 { 1.0 / *d } else { 1.0 });
    let boxed: Box<dyn Any> = Box::new(diag);
    match boxed.downcast::<Vec<S>>() {
        // `S` is `f64`: same type, no copy.
        Ok(same) => Ok(*same),
        Err(other) => {
            let diag: Vec<f64> = match other.downcast::<Vec<f64>>() {
                Ok(v) => *v,
                // Not reachable: the box holds a `Vec<f64>`.
                Err(_) => Vec::new(),
            };
            let mut out: Vec<S> = tracker.try_zeroed(diag.len())?;
            out.par_iter_mut()
                .zip(diag.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(o, d)| *o = S::from_f64(*d));
            tracker.release((diag.len() as u64).saturating_mul(8));
            Ok(out)
        }
    }
}

/// Preconditioned CG (Jacobi) as a state machine.
///
/// Full vectors: `x`, `r`, `p`, `ap` and `inv_diag` (5 x size of `S` per dof).
/// There is no `z` and no `rhs`: `r` starts as `rhs` (the start value `x0`
/// is zero) and `rz = sum(r * r * inv_diag)`. Every sum adds in `f64`.
/// [`restart_from`](Self::restart_from) sets a start value that is not zero.
#[derive(Debug)]
pub struct CgState<S: Scalar> {
    pub x: Vec<S>,
    pub r: Vec<S>,
    pub p: Vec<S>,
    pub ap: Vec<S>,
    pub inv_diag: Vec<S>,
    rz: f64,
    rhs_norm: f64,
    residual: f64,
    tol: f64,
    max_iter: usize,
    iterations: usize,
    done: Option<CgStatus>,
    /// Best residual seen below the stagnation level (single precision only).
    best: f64,
    /// Iterations since `best` fell by 1 percent (single precision only).
    since_best: usize,
    /// Why a single precision run stopped early.
    reason: Option<&'static str>,
}

impl<S: Scalar> CgState<S> {
    /// Make the state. `rhs` becomes `r` (no copy). `diag` is turned into
    /// `inv_diag`. `x`, `p` and `ap` come from `tracker`.
    pub fn new(
        tracker: &mut MemTracker,
        rhs: Vec<S>,
        diag: Vec<f64>,
        tol: f64,
        max_iter: usize,
    ) -> Result<Self, CoreError> {
        let n = rhs.len();
        // Do this first: for `f32` the `f64` diagonal is freed before the
        // other vectors are made.
        let inv_diag: Vec<S> = jacobi_inverse(tracker, diag)?;
        let x: Vec<S> = tracker.try_zeroed(n)?;
        let mut p: Vec<S> = tracker.try_zeroed(n)?;
        let ap: Vec<S> = tracker.try_zeroed(n)?;
        let r = rhs;

        p.par_iter_mut()
            .zip(r.par_iter())
            .zip(inv_diag.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|((pi, ri), di)| *pi = S::from_f64(ri.to_f64() * di.to_f64()));
        let rz = weighted_dot(&r, &inv_diag);
        let r_norm = dot_s(&r, &r).sqrt();
        let rhs_norm = r_norm.max(1e-300);
        let residual = r_norm / rhs_norm;

        let done = if residual <= tol {
            Some(CgStatus::Converged)
        } else if max_iter == 0 {
            Some(CgStatus::MaxIter)
        } else {
            None
        };
        Ok(Self {
            x,
            r,
            p,
            ap,
            inv_diag,
            rz,
            rhs_norm,
            residual,
            tol,
            max_iter,
            iterations: 0,
            done,
            best: f64::INFINITY,
            since_best: 0,
            reason: None,
        })
    }

    /// Start from `x0` and not from zero. Call it before the first `step`.
    /// `r` becomes `rhs - N x0`. `ap` is the work vector for `N x0`, so no
    /// new vector is made. The residual stays relative to `|rhs|`.
    /// `op(v, out)` must write `out = N v`.
    pub fn restart_from(&mut self, x0: &[f64], op: &mut dyn FnMut(&[S], &mut [S])) {
        debug_assert_eq!(self.iterations, 0);
        self.x
            .par_iter_mut()
            .zip(x0.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(xi, v)| *xi = S::from_f64(*v));
        op(&self.x, &mut self.ap); // ap = N x0
        self.r
            .par_iter_mut()
            .zip(self.ap.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(ri, a)| *ri = S::from_f64(ri.to_f64() - a.to_f64()));
        self.p
            .par_iter_mut()
            .zip(self.r.par_iter())
            .zip(self.inv_diag.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|((pi, ri), di)| *pi = S::from_f64(ri.to_f64() * di.to_f64()));
        self.rz = weighted_dot(&self.r, &self.inv_diag);
        self.residual = dot_s(&self.r, &self.r).sqrt() / self.rhs_norm;
        self.done = if self.residual <= self.tol {
            Some(CgStatus::Converged)
        } else if self.max_iter == 0 {
            Some(CgStatus::MaxIter)
        } else {
            None
        };
    }

    /// Number of CG iterations done.
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Relative residual `|r| / |rhs|`.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Text that tells why a single precision run stopped early, if it did.
    pub fn stop_reason(&self) -> Option<&'static str> {
        self.reason
    }

    /// Free `r`, `p`, `ap` and `inv_diag`. Call it only when the run is done.
    /// `tracker` gets the bytes back.
    pub fn release_work(&mut self, tracker: &mut MemTracker) {
        let item = std::mem::size_of::<S>() as u64;
        let count = self.r.len() + self.p.len() + self.ap.len() + self.inv_diag.len();
        self.r = Vec::new();
        self.p = Vec::new();
        self.ap = Vec::new();
        self.inv_diag = Vec::new();
        tracker.release((count as u64).saturating_mul(item));
    }

    /// Run at most `iters` iterations. `op(v, out)` must write `out = N v`.
    pub fn step(&mut self, op: &mut dyn FnMut(&[S], &mut [S]), iters: usize) -> CgStatus {
        if let Some(status) = self.done {
            return status;
        }
        let single = is_single::<S>();
        let mut left = iters;
        while left > 0 {
            if self.iterations >= self.max_iter {
                self.done = Some(CgStatus::MaxIter);
                return CgStatus::MaxIter;
            }
            op(&self.p, &mut self.ap); // ap = N p
            let p_ap = dot_s(&self.p, &self.ap);
            if single && !p_ap.is_finite() {
                self.reason = Some("the residual became NaN");
                self.done = Some(CgStatus::Breakdown);
                return CgStatus::Breakdown;
            }
            if p_ap.abs() < 1e-300 {
                self.done = Some(CgStatus::Breakdown);
                return CgStatus::Breakdown;
            }
            let alpha = self.rz / p_ap;
            self.x
                .par_iter_mut()
                .zip(self.p.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(xi, pi)| *xi = S::from_f64(xi.to_f64() + alpha * pi.to_f64()));
            self.r
                .par_iter_mut()
                .zip(self.ap.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|(ri, api)| *ri = S::from_f64(ri.to_f64() - alpha * api.to_f64()));
            self.iterations += 1;
            left -= 1;
            self.residual = dot_s(&self.r, &self.r).sqrt() / self.rhs_norm;
            if self.residual <= self.tol {
                self.done = Some(CgStatus::Converged);
                return CgStatus::Converged;
            }
            if single {
                if !self.residual.is_finite() {
                    self.reason = Some("the residual became NaN");
                    self.done = Some(CgStatus::Breakdown);
                    return CgStatus::Breakdown;
                }
                if self.residual < F32_STAGNATION_LEVEL {
                    if self.residual < self.best * 0.99 {
                        self.best = self.residual;
                        self.since_best = 0;
                    } else {
                        self.since_best += 1;
                        if self.since_best >= F32_STAGNATION_ITERS {
                            self.reason = Some("the residual does not fall any more");
                            self.done = Some(CgStatus::Breakdown);
                            return CgStatus::Breakdown;
                        }
                    }
                }
            }
            let rz_new = weighted_dot(&self.r, &self.inv_diag);
            let beta = rz_new / self.rz;
            self.p
                .par_iter_mut()
                .zip(self.r.par_iter())
                .zip(self.inv_diag.par_iter())
                .with_min_len(PAR_MIN_LEN)
                .for_each(|((pi, ri), di)| {
                    *pi = S::from_f64(ri.to_f64() * di.to_f64() + beta * pi.to_f64())
                });
            self.rz = rz_new;
        }
        if self.iterations >= self.max_iter {
            self.done = Some(CgStatus::MaxIter);
            return CgStatus::MaxIter;
        }
        CgStatus::Running {
            iteration: self.iterations,
            residual: self.residual,
        }
    }
}

/// `sum(r * r * w)`, added in `f64`.
#[inline]
pub(crate) fn weighted_dot<S: Scalar>(r: &[S], w: &[S]) -> f64 {
    r.par_iter()
        .zip(w.par_iter())
        .with_min_len(PAR_MIN_LEN)
        .map(|(ri, wi)| ri.to_f64() * ri.to_f64() * wi.to_f64())
        .sum::<f64>()
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
