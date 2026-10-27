//! ADMM solver state machine for inequality bounds.
//!
//! The method is the one of `FiniteDifferenceInterpolator::solve_admm`.
//! Outer loop: `z = clip(Q x + u, lower, upper)`, then `u += Q x - z`.
//! Inner loop: Jacobi CG on `(N + rho Q^T Q) x = rhs + rho Q^T (z - u)`.
//!
//! Full-size vectors: `x`, `r`, `p`, `ap` and `inv_diag`.
//! ADMM adds no full-size vector to plain CG. The data right-hand side is
//! computed again from the constraint rows when `r` is rebuilt.
//! `r` is rebuilt at the start of each outer iteration. No dense matrix is made.
//! Vectors with one entry per inequality row: `z`, `u`, `lower`, `upper`.
//! The row order is: value-bound rows, then pair-bound rows.

use crate::PAR_MIN_LEN;
use crate::cg::{dot_s, is_single, jacobi_inverse, weighted_dot};
use crate::error::{BoundKind, CoreError, ViolatedConstraint};
use crate::interpolator::FiniteDifferenceInterpolator;
use crate::memory::MemTracker;
use crate::par::*;
use crate::scalar::Scalar;
use crate::settings::Settings;
use std::any::Any;

/// Number of outer iterations without a 1 percent fall of the primal residual.
const STALL_OUTER: usize = 20;
/// Number of violated rows in an `Infeasible` error.
const MAX_WORST: usize = 20;
/// Largest allowed violation of a bound, as a fraction of the field range
/// (max(x) - min(x)). This is SC-004 of spec 001-low-memory-interpolator:
/// every inequality constraint holds with a violation of at most 1e-4 of the
/// field range. The violation is measured in field units.
pub(crate) const VIOLATION_LIMIT: f64 = 1e-4;

/// Result of one `AdmmState::step` call.
#[derive(Debug, Clone, Copy)]
pub(crate) enum AdmmOutcome {
    /// The solve goes on. `check` is true if an outer iteration just ended.
    Running {
        check: bool,
        iteration: usize,
        residual: f64,
        max_violation: f64,
    },
    /// The solve has ended. The value is true if it converged.
    Done(bool),
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Phase {
    /// The next call must start the pre-solve.
    Start,
    /// The plain CG pre-solve of `N x = rhs` runs.
    Pre,
    /// The next call must start an outer iteration.
    Begin,
    /// The inner CG runs.
    Inner,
}

/// ADMM as a state machine.
#[derive(Debug)]
pub(crate) struct AdmmState<S: Scalar> {
    pub x: Vec<S>,
    r: Vec<S>,
    p: Vec<S>,
    ap: Vec<S>,
    inv_diag: Vec<S>,
    z: Vec<f64>,
    u: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    /// Per row: factor from row units to field units (`1 / scale`).
    inv_scale: Vec<f64>,
    /// Tolerance of the pre-solve (`Settings.tol`).
    tol: f64,
    /// Tolerance of the CG run that is active now.
    cur_tol: f64,
    /// Number of value-bound rows (the first rows).
    n_value: usize,
    rho: f64,
    inner_tol: f64,
    inner_max_iter: usize,
    abs_tol: f64,
    rel_tol: f64,
    /// Limit of the total inner CG iterations.
    max_total: usize,
    total_iters: usize,
    outer_iters: usize,
    // Inner CG state.
    phase: Phase,
    rz: f64,
    rhs_norm: f64,
    inner_res: f64,
    inner_iters: usize,
    breakdown: bool,
    // Outer state.
    primal_rel: f64,
    max_violation: f64,
    best_primal: f64,
    stall: usize,
    done: Option<bool>,
    /// True for single precision: NaN and stagnation tests are on.
    single: bool,
}

/// `out = (N + rho Q^T Q) v`.
fn apply<S: Scalar>(
    interp: &FiniteDifferenceInterpolator,
    scratch: &mut [S],
    rho: f64,
    v: &[S],
    out: &mut [S],
) {
    interp.apply_normal_operator_into(v, scratch, out);
    interp.accumulate_inequality_normal(v, rho, out);
}

/// How far `q` is outside `[lower, upper]`.
#[inline]
fn violation(q: f64, lower: f64, upper: f64) -> f64 {
    (lower - q).max(q - upper).max(0.0)
}

impl<S: Scalar> AdmmState<S> {
    /// Make the state. All large vectors come from `tracker`.
    pub(crate) fn new(
        tracker: &mut MemTracker,
        interp: &FiniteDifferenceInterpolator,
        s: &Settings,
    ) -> Result<Self, CoreError> {
        let n = interp.dof();
        let m = interp.n_inequality_constraints();
        // The diagonal is the one of N only. The pre-solve uses it. The
        // inequality part is added when the ADMM loop starts.
        let diag = interp.try_diagonal(tracker)?;
        // Jacobi: a diagonal value of zero or less gives no scaling.
        // For f32 the f64 diagonal is freed before the other vectors are made.
        let inv_diag: Vec<S> = jacobi_inverse(tracker, diag)?;
        let x: Vec<S> = tracker.try_zeroed(n)?;
        let r: Vec<S> = tracker.try_zeroed(n)?;
        let p: Vec<S> = tracker.try_zeroed(n)?;
        let ap: Vec<S> = tracker.try_zeroed(n)?;
        let z: Vec<f64> = tracker.try_zeroed(m)?;
        let u: Vec<f64> = tracker.try_zeroed(m)?;
        let (lower, upper) = interp.inequality_bounds();
        Ok(Self {
            x,
            r,
            p,
            ap,
            inv_diag,
            z,
            u,
            lower,
            upper,
            inv_scale: interp.inequality_inverse_scales(),
            tol: s.tol,
            cur_tol: s.tol,
            n_value: interp.n_value_inequality_constraints(),
            rho: s.rho,
            // The inner solve is at least as tight as the global tolerance.
            inner_tol: s.inner_tol.min(s.tol),
            inner_max_iter: s.inner_max_iter,
            abs_tol: s.admm_abs_tol,
            rel_tol: s.admm_rel_tol,
            max_total: s.max_iter,
            total_iters: 0,
            outer_iters: 0,
            phase: Phase::Start,
            rz: 0.0,
            rhs_norm: 1e-300,
            inner_res: 1.0,
            inner_iters: 0,
            breakdown: false,
            primal_rel: 1.0,
            max_violation: 0.0,
            best_primal: f64::INFINITY,
            stall: 0,
            done: None,
            single: is_single::<S>(),
        })
    }

    /// Free `r`, `p`, `ap` and `inv_diag`. Call it only when the run is done.
    pub(crate) fn release_work(&mut self, tracker: &mut MemTracker) {
        let item = std::mem::size_of::<S>() as u64;
        let count = self.r.len() + self.p.len() + self.ap.len() + self.inv_diag.len();
        self.r = Vec::new();
        self.p = Vec::new();
        self.ap = Vec::new();
        self.inv_diag = Vec::new();
        tracker.release((count as u64).saturating_mul(item));
    }

    /// Total inner CG iterations done.
    pub(crate) fn iterations(&self) -> usize {
        self.total_iters
    }

    /// Relative primal residual of the last outer iteration.
    pub(crate) fn residual(&self) -> f64 {
        self.primal_rel
    }

    /// Largest violation of the bounds at the last outer iteration end.
    pub(crate) fn max_violation(&self) -> f64 {
        self.max_violation
    }

    /// Start an outer iteration: build `r = rhs_admm - A x` and the first `p`.
    fn begin_outer(&mut self, interp: &FiniteDifferenceInterpolator, scratch: &mut [S]) {
        let m = self.z.len();
        // `zu` is small: one entry per inequality row.
        let zu: Vec<f64> = (0..m).map(|i| self.z[i] - self.u[i]).collect();
        // Data right-hand side, computed again into `r`.
        interp.write_data_rhs(&mut self.r);
        interp.accumulate_inequality_transpose(&zu, self.rho, &mut self.r);
        // The inner tolerance is relative to the norm of the ADMM right-hand side.
        self.rhs_norm = dot_s(&self.r, &self.r).sqrt().max(1e-300);
        apply(interp, scratch, self.rho, &self.x, &mut self.ap);
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
        self.inner_res = dot_s(&self.r, &self.r).sqrt() / self.rhs_norm;
        self.inner_iters = 0;
        self.breakdown = false;
        self.cur_tol = self.inner_tol;
        self.phase = Phase::Inner;
    }

    /// Start the pre-solve: plain CG of `N x = rhs` from `x = 0`.
    /// `r` is the data right-hand side (x = 0). No new vector is needed.
    fn begin_pre(&mut self, interp: &FiniteDifferenceInterpolator) {
        interp.write_data_rhs(&mut self.r);
        self.p
            .par_iter_mut()
            .zip(self.r.par_iter())
            .zip(self.inv_diag.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|((pi, ri), di)| *pi = S::from_f64(ri.to_f64() * di.to_f64()));
        self.rz = weighted_dot(&self.r, &self.inv_diag);
        let r_norm = dot_s(&self.r, &self.r).sqrt();
        self.rhs_norm = r_norm.max(1e-300);
        self.inner_res = r_norm / self.rhs_norm;
        self.inner_iters = 0;
        self.breakdown = false;
        self.cur_tol = self.tol;
        self.phase = Phase::Pre;
    }

    /// `max(x) - min(x)`, with a floor so that the limit is never zero.
    fn field_range(&self) -> f64 {
        let (lo, hi) = self
            .x
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                let v = v.to_f64();
                (lo.min(v), hi.max(v))
            });
        (hi - lo).max(1e-12)
    }

    /// Largest violation of the bounds in field units.
    fn field_violation(&self, qx: &[f64]) -> f64 {
        let mut m = 0.0f64;
        for (i, q) in qx.iter().enumerate() {
            m = m.max(violation(*q, self.lower[i], self.upper[i]) * self.inv_scale[i]);
        }
        m
    }

    /// End of the pre-solve. Returns true if the plain solution already
    /// holds all bounds. If not, it makes the ADMM start from this `x`.
    fn end_pre(
        &mut self,
        tracker: &mut MemTracker,
        interp: &FiniteDifferenceInterpolator,
    ) -> Result<bool, CoreError> {
        let qx = interp.apply_inequality_operator(&self.x);
        self.max_violation = self.field_violation(&qx);
        self.primal_rel = self.inner_res;
        if self.max_violation <= VIOLATION_LIMIT * self.field_range() {
            return Ok(true);
        }
        // Warm start: z = clip(Q x), u = 0.
        for (i, q) in qx.iter().enumerate() {
            self.z[i] = q.max(self.lower[i]).min(self.upper[i]);
            self.u[i] = 0.0;
        }
        // Put the inequality part into the Jacobi scaling.
        // An entry that was zero or less was stored as 1; it stays a weak
        // scaling only. This changes the speed, not the result.
        if let Some(d) = (&mut self.inv_diag as &mut dyn Any).downcast_mut::<Vec<f64>>() {
            // Double precision: the change is done in place.
            d.par_iter_mut()
                .with_min_len(PAR_MIN_LEN)
                .for_each(|d| *d = 1.0 / *d);
            interp.add_inequality_diagonal(self.rho, d);
            d.par_iter_mut()
                .with_min_len(PAR_MIN_LEN)
                .for_each(|d| *d = if *d > 0.0 { 1.0 / *d } else { 1.0 });
            return Ok(false);
        }
        // Single precision: `add_inequality_diagonal` needs f64 values.
        // Use one short-lived f64 vector, counted by the tracker.
        let n = self.inv_diag.len();
        let mut tmp: Vec<f64> = tracker.try_zeroed(n)?;
        tmp.par_iter_mut()
            .zip(self.inv_diag.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(t, d)| *t = 1.0 / d.to_f64());
        interp.add_inequality_diagonal(self.rho, &mut tmp);
        self.inv_diag
            .par_iter_mut()
            .zip(tmp.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(d, t)| *d = S::from_f64(if *t > 0.0 { 1.0 / *t } else { 1.0 }));
        tracker.release((n as u64).saturating_mul(8));
        Ok(false)
    }

    /// One inner CG iteration.
    fn cg_iteration(&mut self, interp: &FiniteDifferenceInterpolator, scratch: &mut [S]) {
        if self.phase == Phase::Pre {
            interp.apply_normal_operator_into(&self.p, scratch, &mut self.ap);
        } else {
            apply(interp, scratch, self.rho, &self.p, &mut self.ap);
        }
        let p_ap = dot_s(&self.p, &self.ap);
        if p_ap.abs() < 1e-300 || (self.single && !p_ap.is_finite()) {
            self.breakdown = true;
            return;
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
        self.total_iters += 1;
        self.inner_iters += 1;
        self.inner_res = dot_s(&self.r, &self.r).sqrt() / self.rhs_norm;
        if self.inner_res <= self.cur_tol {
            return;
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

    fn inner_finished(&self) -> bool {
        let own_limit = self.phase != Phase::Pre && self.inner_iters >= self.inner_max_iter;
        self.breakdown
            || (self.single && !self.inner_res.is_finite())
            || self.inner_res <= self.cur_tol
            || own_limit
            || self.total_iters >= self.max_total
    }

    /// The `z` and `u` update after the inner solve. Returns true if converged.
    fn end_outer(&mut self, interp: &FiniteDifferenceInterpolator) -> Result<bool, CoreError> {
        let m = self.z.len();
        let qx = interp.apply_inequality_operator(&self.x);
        let mut primal2 = 0.0f64;
        let mut dual2 = 0.0f64;
        let mut qx2 = 0.0f64;
        let mut z2 = 0.0f64;
        let mut u2 = 0.0f64;
        let mut max_violation = 0.0f64;
        for i in 0..m {
            let v = qx[i] + self.u[i];
            let z_new = v.max(self.lower[i]).min(self.upper[i]);
            dual2 += (z_new - self.z[i]).powi(2);
            self.z[i] = z_new;
            self.u[i] += qx[i] - z_new;
            primal2 += (qx[i] - z_new).powi(2);
            qx2 += qx[i] * qx[i];
            z2 += z_new * z_new;
            u2 += self.u[i] * self.u[i];
            // Field units: undo the row scaling.
            max_violation = max_violation
                .max(violation(qx[i], self.lower[i], self.upper[i]) * self.inv_scale[i]);
        }
        self.max_violation = max_violation;
        let limit = VIOLATION_LIMIT * self.field_range();
        self.outer_iters += 1;

        let primal_norm = primal2.sqrt();
        let dual_norm = self.rho * dual2.sqrt();
        let qx_norm = qx2.sqrt();
        let z_norm = z2.sqrt();
        let u_norm = self.rho * u2.sqrt();
        let sqrt_m = (m as f64).sqrt();
        let eps_pri = sqrt_m * self.abs_tol + self.rel_tol * qx_norm.max(z_norm);
        let eps_dual = sqrt_m * self.abs_tol + self.rel_tol * u_norm;
        self.primal_rel = primal_norm / qx_norm.max(z_norm).max(1e-300);

        // Converged only if the bounds also hold to SC-004.
        if primal_norm <= eps_pri && dual_norm <= eps_dual && max_violation <= limit {
            return Ok(true);
        }

        // Stall test: the primal residual must fall by 1 percent in 20 outer iterations.
        if primal_norm < self.best_primal * 0.99 {
            self.best_primal = primal_norm;
            self.stall = 0;
        } else {
            self.stall += 1;
        }
        if self.stall >= STALL_OUTER && max_violation > limit {
            return Err(CoreError::Infeasible {
                worst: self.worst_violations(&qx),
            });
        }
        Ok(false)
    }

    /// The largest violations, largest first. At most 20 entries.
    fn worst_violations(&self, qx: &[f64]) -> Vec<ViolatedConstraint> {
        let mut all: Vec<ViolatedConstraint> = Vec::new();
        for (i, q) in qx.iter().enumerate() {
            let v = violation(*q, self.lower[i], self.upper[i]) * self.inv_scale[i];
            if v > 0.0 {
                let (kind, index) = if i < self.n_value {
                    (BoundKind::Value, i)
                } else {
                    (BoundKind::Pair, i - self.n_value)
                };
                all.push(ViolatedConstraint {
                    index,
                    kind,
                    violation: v,
                });
            }
        }
        all.sort_by(|a, b| {
            b.violation
                .partial_cmp(&a.violation)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        all.truncate(MAX_WORST);
        all
    }

    /// Run at most `max_iters` inner CG iterations. The call also returns at
    /// the end of each outer iteration.
    pub(crate) fn step(
        &mut self,
        tracker: &mut MemTracker,
        interp: &FiniteDifferenceInterpolator,
        scratch: &mut [S],
        max_iters: usize,
    ) -> Result<AdmmOutcome, CoreError> {
        if let Some(c) = self.done {
            return Ok(AdmmOutcome::Done(c));
        }
        if self.phase == Phase::Start {
            self.begin_pre(interp);
        }
        if self.phase == Phase::Pre {
            let mut left = max_iters;
            while left > 0 && !self.inner_finished() {
                self.cg_iteration(interp, scratch);
                left -= 1;
            }
            if !self.inner_finished() {
                return Ok(AdmmOutcome::Running {
                    check: false,
                    iteration: self.total_iters,
                    residual: self.inner_res,
                    max_violation: self.max_violation,
                });
            }
            if self.inner_res > self.tol || (self.single && !self.inner_res.is_finite()) {
                // The pre-solve did not reach its tolerance.
                self.primal_rel = self.inner_res;
                self.done = Some(false);
                return Ok(AdmmOutcome::Done(false));
            }
            if self.end_pre(tracker, interp)? {
                self.done = Some(true);
                return Ok(AdmmOutcome::Done(true));
            }
            self.phase = Phase::Begin;
            return Ok(AdmmOutcome::Running {
                check: true,
                iteration: self.total_iters,
                residual: self.primal_rel,
                max_violation: self.max_violation,
            });
        }
        if self.phase == Phase::Begin {
            if self.total_iters >= self.max_total || self.outer_iters >= self.max_total {
                self.done = Some(false);
                return Ok(AdmmOutcome::Done(false));
            }
            self.begin_outer(interp, scratch);
        }
        let mut left = max_iters;
        while left > 0 && !self.inner_finished() {
            self.cg_iteration(interp, scratch);
            left -= 1;
        }
        if !self.inner_finished() {
            return Ok(AdmmOutcome::Running {
                check: false,
                iteration: self.total_iters,
                residual: self.primal_rel,
                max_violation: self.max_violation,
            });
        }
        self.phase = Phase::Begin;
        if self.single && !self.inner_res.is_finite() {
            // Single precision: a NaN residual ends the solve.
            self.done = Some(false);
            return Ok(AdmmOutcome::Done(false));
        }
        if self.end_outer(interp)? {
            self.done = Some(true);
            return Ok(AdmmOutcome::Done(true));
        }
        if self.total_iters >= self.max_total {
            self.done = Some(false);
            return Ok(AdmmOutcome::Done(false));
        }
        Ok(AdmmOutcome::Running {
            check: true,
            iteration: self.total_iters,
            residual: self.primal_rel,
            max_violation: self.max_violation,
        })
    }
}
