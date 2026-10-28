//! Solve session: a state machine that gives control back often.

use crate::admm::{AdmmOutcome, AdmmState};
use crate::cg::{CgState, CgStatus};
use crate::error::CoreError;
use crate::field::{Field, FieldData};
use crate::interpolator::FiniteDifferenceInterpolator;
use crate::memory::MemTracker;
use crate::problem::Problem;
use crate::report::{Progress, SolveReport, Stage, Step};
use crate::scalar::Scalar;
use crate::settings::{Precision, Settings};

/// The solver of a session: plain CG, or ADMM if the problem has bounds.
#[derive(Debug)]
enum Solver<S: Scalar> {
    Cg(CgState<S>),
    Admm(AdmmState<S>),
}

/// Result of one step of a `State`.
enum Out {
    Running {
        check: bool,
        iteration: usize,
        residual: f64,
        max_violation: f64,
    },
    /// The solve has ended. The value is true if the solver says it converged.
    Done(bool),
}

/// Solver and scratch buffer, with the vector type `S` (`f64` or `f32`).
#[derive(Debug)]
struct State<S: Scalar> {
    solver: Solver<S>,
    /// Scratch of the normal operator. It is allocated once.
    scratch: Vec<S>,
}

impl<S: Scalar> State<S> {
    /// Allocate the solver vectors and the scratch buffer from `tracker`.
    /// With `x0` the solver starts from `x0`, not from zero.
    fn new(
        interp: &FiniteDifferenceInterpolator,
        tracker: &mut MemTracker,
        s: &Settings,
        x0: Option<&[f64]>,
    ) -> Result<Self, CoreError> {
        // The `f32` scratch has no border part if the problem needs none.
        // The `f64` scratch is not changed.
        let scratch_len = if std::mem::size_of::<S>() == 4 {
            interp.scratch_len_compact()
        } else {
            interp.scratch_len()
        };
        let solver = if interp.n_inequality_constraints() > 0 {
            Solver::Admm(AdmmState::new(tracker, interp, s)?)
        } else {
            // `rhs` becomes `r` in the CG state. It is not copied.
            let mut rhs: Vec<S> = tracker.try_zeroed(interp.dof())?;
            interp.write_data_rhs(&mut rhs);
            let diag = interp.try_diagonal(tracker)?;
            Solver::Cg(CgState::new(tracker, rhs, diag, s.tol, s.max_iter)?)
        };
        let mut scratch: Vec<S> = tracker.try_zeroed(scratch_len)?;
        let mut solver = solver;
        if let Some(x0) = x0 {
            match &mut solver {
                Solver::Cg(cg) => {
                    let mut op = |v: &[S], out: &mut [S]| {
                        interp.apply_normal_operator_into(v, &mut scratch, out);
                    };
                    cg.restart_from(x0, &mut op);
                }
                Solver::Admm(admm) => admm.restart_from(x0),
            }
        }
        Ok(Self { solver, scratch })
    }

    fn iterations_residual_violation(&self) -> (usize, f64, f64) {
        match &self.solver {
            Solver::Cg(cg) => (cg.iterations(), cg.residual(), 0.0),
            Solver::Admm(a) => (a.iterations(), a.residual(), a.max_violation()),
        }
    }

    fn x(&self) -> &[S] {
        match &self.solver {
            Solver::Cg(cg) => &cg.x,
            Solver::Admm(a) => &a.x,
        }
    }

    fn x_mut(&mut self) -> &mut Vec<S> {
        match &mut self.solver {
            Solver::Cg(cg) => &mut cg.x,
            Solver::Admm(a) => &mut a.x,
        }
    }

    fn into_x(self) -> Vec<S> {
        match self.solver {
            Solver::Cg(cg) => cg.x,
            Solver::Admm(a) => a.x,
        }
    }

    fn is_cg(&self) -> bool {
        matches!(self.solver, Solver::Cg(_))
    }

    /// Text that tells why a single precision CG run stopped early.
    fn stop_reason(&self) -> Option<&'static str> {
        match &self.solver {
            Solver::Cg(cg) => cg.stop_reason(),
            Solver::Admm(_) => None,
        }
    }

    /// Free all work vectors and the scratch buffer. Only `x` stays.
    /// Call it only when the run is done.
    fn release_work(&mut self, tracker: &mut MemTracker) {
        match &mut self.solver {
            Solver::Cg(cg) => cg.release_work(tracker),
            Solver::Admm(a) => a.release_work(tracker),
        }
        let bytes = (self.scratch.len() as u64).saturating_mul(std::mem::size_of::<S>() as u64);
        self.scratch = Vec::new();
        tracker.release(bytes);
    }

    /// Run at most `max_iters` iterations.
    fn step(
        &mut self,
        interp: &FiniteDifferenceInterpolator,
        tracker: &mut MemTracker,
        max_iters: usize,
    ) -> Result<Out, CoreError> {
        let scratch = &mut self.scratch;
        match &mut self.solver {
            Solver::Cg(cg) => {
                let mut op = |v: &[S], out: &mut [S]| {
                    interp.apply_normal_operator_into(v, scratch, out);
                };
                Ok(match cg.step(&mut op, max_iters) {
                    CgStatus::Running {
                        iteration,
                        residual,
                    } => Out::Running {
                        check: false,
                        iteration,
                        residual,
                        max_violation: 0.0,
                    },
                    CgStatus::Converged => Out::Done(true),
                    CgStatus::MaxIter | CgStatus::Breakdown => Out::Done(false),
                })
            }
            Solver::Admm(admm) => Ok(match admm.step(tracker, interp, scratch, max_iters)? {
                AdmmOutcome::Running {
                    check,
                    iteration,
                    residual,
                    max_violation,
                } => Out::Running {
                    check,
                    iteration,
                    residual,
                    max_violation,
                },
                AdmmOutcome::Done(converged) => Out::Done(converged),
            }),
        }
    }
}

/// The state of a session, in the precision of `Settings.precision`.
#[derive(Debug)]
enum Inner {
    F64(State<f64>),
    F32(State<f32>),
}

/// Copy a slice into a new vector. A failed allocation is an error.
fn copy_vec<S: Scalar>(src: &[S]) -> Result<Vec<S>, CoreError> {
    let mut x: Vec<S> = Vec::new();
    x.try_reserve_exact(src.len())
        .map_err(|_| CoreError::OutOfMemory {
            requested_bytes: (src.len() as u64).saturating_mul(std::mem::size_of::<S>() as u64),
        })?;
    x.extend_from_slice(src);
    Ok(x)
}

/// A running solve. The session owns the problem. Drop it to cancel.
#[derive(Debug)]
pub struct SolveSession {
    problem: Problem,
    inner: Inner,
    tracker: MemTracker,
    stage: Stage,
    /// Set when the solve has ended: true if it converged.
    done: Option<bool>,
    /// Tolerance of the solve.
    tol: f64,
    /// Single precision CG: the residual from `x`, computed again in `f64`.
    final_residual: Option<f64>,
    /// Warnings that the solver adds to the report.
    notes: Vec<String>,
    /// Set when `finalize` fails. The session is then broken: every later
    /// `step` returns this error. The state is not valid any more.
    failed: Option<CoreError>,
    /// Test only: make the next `finalize` fail once with this error.
    #[cfg(test)]
    inject_finalize_error: Option<CoreError>,
}

impl SolveSession {
    /// Set up the solve. The stage is `Setup` until the first `step`.
    pub(crate) fn new(mut problem: Problem, s: &Settings) -> Result<Self, CoreError> {
        let mut tracker = MemTracker::new(s.budget_bytes);
        problem.interpolator_mut().prepare(&mut tracker)?;
        // The solver copies the start value. Then it is freed.
        let x0 = problem.start.take();
        let x0 = x0.as_deref();
        let interp = problem.interpolator();
        let inner = match s.precision {
            Precision::F64 => Inner::F64(State::new(interp, &mut tracker, s, x0)?),
            Precision::F32 => Inner::F32(State::new(interp, &mut tracker, s, x0)?),
        };
        Ok(Self {
            problem,
            inner,
            tracker,
            stage: Stage::Setup,
            done: None,
            tol: s.tol,
            final_residual: None,
            notes: Vec::new(),
            failed: None,
            #[cfg(test)]
            inject_finalize_error: None,
        })
    }

    fn report(&self, converged: bool) -> SolveReport {
        let (iterations, residual, max_violation) = match &self.inner {
            Inner::F64(st) => st.iterations_residual_violation(),
            Inner::F32(st) => st.iterations_residual_violation(),
        };
        let mut warnings = self.problem.warnings.clone();
        warnings.extend(self.notes.iter().cloned());
        SolveReport {
            converged,
            iterations,
            residual: self.final_residual.unwrap_or(residual),
            max_violation,
            peak_bytes: self.tracker.peak(),
            warnings,
        }
    }

    /// End of a single precision solve. For plain CG, free the work vectors
    /// and compute the residual of the final `x` again in `f64`. Add a
    /// warning if the solve did not converge. Returns the final `converged`.
    /// A double precision solve is not changed.
    fn finalize(&mut self, converged: bool) -> Result<bool, CoreError> {
        let Inner::F32(st) = &mut self.inner else {
            return Ok(converged);
        };
        let reason = st.stop_reason();
        let mut ok = converged;
        let mut true_residual: Option<f64> = None;
        // The work vectors are not needed any more. Free them first, so the
        // peak of the check stays low.
        st.release_work(&mut self.tracker);
        #[cfg(test)]
        if let Some(e) = self.inject_finalize_error.take() {
            return Err(e);
        }
        if st.is_cg() {
            let n = st.x().len();
            // One `f64` copy of `x` is needed for the check. If it does not
            // fit, keep the residual of the iteration and say so.
            match self.tracker.try_zeroed::<f64>(n) {
                Ok(mut x64) => {
                    for (d, s) in x64.iter_mut().zip(st.x().iter()) {
                        *d = *s as f64;
                    }
                    // Free the `f32` copy during the check, then restore it.
                    // The `f32` to `f64` to `f32` trip is exact.
                    // The `f32` buffer cannot be reserved before the check:
                    // that would raise the peak by 4 bytes per dof. If the
                    // restore fails, `step` stores the error and the session
                    // stays failed.
                    let x32 = std::mem::take(st.x_mut());
                    let bytes32 = (x32.len() as u64).saturating_mul(4);
                    drop(x32);
                    self.tracker.release(bytes32);
                    let res = self.problem.interpolator().residual_f64_compact(&x64);
                    let mut back: Vec<f32> = self.tracker.try_zeroed(n)?;
                    for (d, s) in back.iter_mut().zip(x64.iter()) {
                        *d = *s as f32;
                    }
                    self.tracker.release((n as u64).saturating_mul(8));
                    drop(x64);
                    *st.x_mut() = back;
                    true_residual = Some(res);
                }
                Err(_) => {
                    self.notes.push(
                        "Single precision: there was no memory to check the residual in \
                         double precision. The residual is the one of the iteration."
                            .into(),
                    );
                }
            }
        }
        if let Some(res) = true_residual {
            self.final_residual = Some(res);
            ok = res <= self.tol;
        }
        if !ok {
            let residual = self
                .final_residual
                .unwrap_or_else(|| st.iterations_residual_violation().1);
            let mut text = format!(
                "Single precision (f32) did not reach the tolerance {:e}. \
                 The residual is {:e}. Use double precision or a larger tolerance.",
                self.tol, residual
            );
            if let Some(r) = reason {
                text.push_str(&format!(" Stop reason: {r}."));
            }
            self.notes.push(text);
        }
        Ok(ok)
    }

    /// Run at most `max_iters` iterations.
    pub fn step(&mut self, max_iters: usize) -> Result<Step, CoreError> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        if let Some(converged) = self.done {
            return Ok(Step::Done(self.report(converged)));
        }
        self.stage = Stage::Solve;
        let interp = self.problem.interpolator();
        let out = match &mut self.inner {
            Inner::F64(st) => st.step(interp, &mut self.tracker, max_iters)?,
            Inner::F32(st) => st.step(interp, &mut self.tracker, max_iters)?,
        };
        match out {
            Out::Running {
                check,
                iteration,
                residual,
                max_violation,
            } => {
                if check {
                    self.stage = Stage::CheckConstraints;
                }
                Ok(Step::Running(Progress {
                    stage: self.stage.clone(),
                    iteration,
                    residual,
                    max_violation,
                }))
            }
            Out::Done(converged) => {
                let converged = match self.finalize(converged) {
                    Ok(c) => c,
                    Err(e) => {
                        self.failed = Some(e.clone());
                        return Err(e);
                    }
                };
                self.done = Some(converged);
                Ok(Step::Done(self.report(converged)))
            }
        }
    }

    /// Return the current field. `converged` is false until `Done`.
    pub fn best_so_far(&self) -> Result<(Field, SolveReport), CoreError> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        let data = match &self.inner {
            Inner::F64(st) => FieldData::F64(copy_vec(st.x())?),
            Inner::F32(st) => FieldData::F32(copy_vec(st.x())?),
        };
        let field = Field {
            grid: self.problem.interpolator().grid,
            data,
        };
        Ok((field, self.report(self.done == Some(true))))
    }

    /// Move the field out.
    pub fn finish(self) -> Result<(Field, SolveReport), CoreError> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        let report = self.report(self.done == Some(true));
        let grid = self.problem.interpolator().grid;
        let data = match self.inner {
            Inner::F64(st) => FieldData::F64(st.into_x()),
            Inner::F32(st) => FieldData::F32(st.into_x()),
        };
        Ok((Field { grid, data }, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f32_session() -> SolveSession {
        let mut p = Problem::new([0.0; 3], [5; 3], [1.0; 3]).unwrap();
        let pts = [1.0, 1.0, 1.0, 3.0, 2.0, 2.0, 2.0, 3.0, 1.0];
        p.add_value(&pts, &[1.0, 2.0, 3.0], &[1.0; 3]).unwrap();
        let s = Settings {
            precision: Precision::F32,
            budget_bytes: Some(1 << 24),
            ..Settings::default()
        };
        p.start(s).unwrap()
    }

    #[test]
    fn failed_finalize_is_final() {
        let mut session = f32_session();
        let err = CoreError::OutOfMemory { requested_bytes: 1 };
        session.inject_finalize_error = Some(err.clone());
        let mut first = None;
        for _ in 0..100_000 {
            match session.step(1 << 20) {
                Ok(Step::Running(_)) => continue,
                other => {
                    first = Some(other);
                    break;
                }
            }
        }
        assert_eq!(first.unwrap().err(), Some(err.clone()));
        // The second call returns the same error. It is not a done step.
        assert_eq!(session.step(1 << 20).err(), Some(err));
        assert!(session.done.is_none());
        // No field of zero length leaves the failed session.
        assert!(session.best_so_far().is_err());
        assert!(session.finish().is_err());
    }

    #[test]
    fn f32_check_keeps_the_field_whole() {
        let mut session = f32_session();
        loop {
            if let Step::Done(_) = session.step(1 << 20).unwrap() {
                break;
            }
        }
        let (field, _) = session.finish().unwrap();
        assert_eq!(field.values_f64().len(), 125);
    }
}
