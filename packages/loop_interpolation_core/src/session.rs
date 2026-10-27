//! Solve session: a state machine that gives control back often.

use crate::admm::{AdmmOutcome, AdmmState};
use crate::cg::{CgState, CgStatus};
use crate::error::CoreError;
use crate::field::{Field, FieldData};
use crate::memory::MemTracker;
use crate::problem::Problem;
use crate::report::{Progress, SolveReport, Stage, Step};
use crate::settings::Settings;

/// The solver of a session: plain CG, or ADMM if the problem has bounds.
#[derive(Debug)]
enum Solver {
    Cg(CgState),
    Admm(AdmmState),
}

/// A running solve. The session owns the problem. Drop it to cancel.
#[derive(Debug)]
pub struct SolveSession {
    problem: Problem,
    solver: Solver,
    /// Scratch of the normal operator. It is allocated once.
    scratch: Vec<f64>,
    tracker: MemTracker,
    stage: Stage,
    /// Set when the solve has ended: true if it converged.
    done: Option<bool>,
}

impl SolveSession {
    /// Set up the solve. The stage is `Setup` until the first `step`.
    pub(crate) fn new(mut problem: Problem, s: &Settings) -> Result<Self, CoreError> {
        let mut tracker = MemTracker::new(s.budget_bytes);
        problem.interpolator_mut().prepare(&mut tracker)?;
        let scratch_len = problem.interpolator().scratch_len();
        let solver = if problem.interpolator().n_inequality_constraints() > 0 {
            Solver::Admm(AdmmState::new(&mut tracker, problem.interpolator(), s)?)
        } else {
            let interp = problem.interpolator();
            // `rhs` becomes `r` in the CG state. It is not copied.
            let rhs = interp.try_rhs(&mut tracker)?;
            let diag = interp.try_diagonal(&mut tracker)?;
            Solver::Cg(CgState::new(&mut tracker, rhs, diag, s.tol, s.max_iter)?)
        };
        let scratch: Vec<f64> = tracker.try_zeroed(scratch_len)?;
        Ok(Self {
            problem,
            solver,
            scratch,
            tracker,
            stage: Stage::Setup,
            done: None,
        })
    }

    fn report(&self, converged: bool) -> SolveReport {
        let (iterations, residual, max_violation) = match &self.solver {
            Solver::Cg(cg) => (cg.iterations(), cg.residual(), 0.0),
            Solver::Admm(a) => (a.iterations(), a.residual(), a.max_violation()),
        };
        SolveReport {
            converged,
            iterations,
            residual,
            max_violation,
            peak_bytes: self.tracker.peak(),
            warnings: self.problem.warnings.clone(),
        }
    }

    /// Run at most `max_iters` iterations.
    pub fn step(&mut self, max_iters: usize) -> Result<Step, CoreError> {
        if let Some(converged) = self.done {
            return Ok(Step::Done(self.report(converged)));
        }
        self.stage = Stage::Solve;
        let interp = self.problem.interpolator();
        let scratch = &mut self.scratch;
        match &mut self.solver {
            Solver::Cg(cg) => {
                let mut op = |v: &[f64], out: &mut [f64]| {
                    interp.apply_normal_operator_into(v, scratch, out);
                };
                let status = cg.step(&mut op, max_iters);
                match status {
                    CgStatus::Running {
                        iteration,
                        residual,
                    } => Ok(Step::Running(Progress {
                        stage: self.stage.clone(),
                        iteration,
                        residual,
                        max_violation: 0.0,
                    })),
                    CgStatus::Converged => {
                        self.done = Some(true);
                        Ok(Step::Done(self.report(true)))
                    }
                    CgStatus::MaxIter | CgStatus::Breakdown => {
                        self.done = Some(false);
                        Ok(Step::Done(self.report(false)))
                    }
                }
            }
            Solver::Admm(admm) => match admm.step(interp, scratch, max_iters)? {
                AdmmOutcome::Running {
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
                AdmmOutcome::Done(converged) => {
                    self.done = Some(converged);
                    Ok(Step::Done(self.report(converged)))
                }
            },
        }
    }

    /// Return the current field. `converged` is false until `Done`.
    pub fn best_so_far(&self) -> Result<(Field, SolveReport), CoreError> {
        let src: &[f64] = match &self.solver {
            Solver::Cg(cg) => &cg.x,
            Solver::Admm(a) => &a.x,
        };
        let mut x: Vec<f64> = Vec::new();
        x.try_reserve_exact(src.len())
            .map_err(|_| CoreError::OutOfMemory {
                requested_bytes: (src.len() as u64).saturating_mul(8),
            })?;
        x.extend_from_slice(src);
        let field = Field {
            grid: self.problem.interpolator().grid,
            data: FieldData::F64(x),
        };
        Ok((field, self.report(self.done == Some(true))))
    }

    /// Move the field out.
    pub fn finish(self) -> Result<(Field, SolveReport), CoreError> {
        let report = self.report(self.done == Some(true));
        let grid = self.problem.interpolator().grid;
        let field = Field {
            grid,
            data: FieldData::F64(match self.solver {
                Solver::Cg(cg) => cg.x,
                Solver::Admm(a) => a.x,
            }),
        };
        Ok((field, report))
    }
}
