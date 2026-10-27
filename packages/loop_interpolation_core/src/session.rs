//! Solve session: a state machine that gives control back often.

use crate::error::CoreError;
use crate::field::Field;
use crate::problem::Problem;
use crate::report::{SolveReport, Step};

/// A running solve. The session owns the problem. Drop it to cancel.
#[derive(Debug, Clone)]
pub struct SolveSession {
    #[allow(dead_code)]
    pub(crate) problem: Problem,
}

fn ni<T>() -> Result<T, CoreError> {
    Err(CoreError::Unsupported {
        setting: "not implemented".into(),
    })
}

impl SolveSession {
    /// Run at most `max_iters` iterations.
    pub fn step(&mut self, _max_iters: usize) -> Result<Step, CoreError> {
        ni()
    }

    /// Return the current field. `converged` is false until `Done`.
    pub fn best_so_far(&self) -> Result<(Field, SolveReport), CoreError> {
        ni()
    }

    /// Move the field out.
    pub fn finish(self) -> Result<(Field, SolveReport), CoreError> {
        ni()
    }
}
