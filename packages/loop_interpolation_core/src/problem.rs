//! The problem: grid, constraints and weights.

use crate::error::CoreError;
use crate::field::Field;
use crate::grid::StructuredGrid3D;
use crate::interpolator::{BorderRegularisationWeights, RegularisationWeights};
use crate::report::{MemoryEstimate, SolveReport};
use crate::session::SolveSession;
use crate::settings::Settings;

fn ni<T>() -> Result<T, CoreError> {
    Err(CoreError::Unsupported {
        setting: "not implemented".into(),
    })
}

/// A flat-input interpolation problem. Arrays are flat: points have 3n values.
#[derive(Debug, Clone)]
pub struct Problem {
    #[allow(dead_code)]
    pub(crate) grid: StructuredGrid3D,
    pub(crate) warnings: Vec<String>,
}

impl Problem {
    pub fn new(origin: [f64; 3], counts: [usize; 3], step: [f64; 3]) -> Result<Problem, CoreError> {
        let _ = (origin, counts, step);
        ni()
    }

    pub fn set_regularisation(&mut self, _w: RegularisationWeights) {
        unimplemented!()
    }

    pub fn set_border_regularisation(&mut self, _w: BorderRegularisationWeights) {
        unimplemented!()
    }

    pub fn add_value(
        &mut self,
        _points: &[f64],
        _values: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    /// Gradient constraint, direction only.
    pub fn add_gradient(
        &mut self,
        _points: &[f64],
        _vectors: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    /// Normal constraint, with length.
    pub fn add_normal(
        &mut self,
        _points: &[f64],
        _vectors: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn add_tangent(
        &mut self,
        _points: &[f64],
        _vectors: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn add_interface(
        &mut self,
        _points: &[f64],
        _ids: &[f64],
        _weight: f64,
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn add_directional_regularisation(
        &mut self,
        _directions: &[f64],
        _weight: f64,
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn add_value_bound(
        &mut self,
        _points: &[f64],
        _lower: &[f64],
        _upper: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn add_pair_bound(
        &mut self,
        _points_a: &[f64],
        _points_b: &[f64],
        _lower: &[f64],
        _upper: &[f64],
        _weights: &[f64],
    ) -> Result<(), CoreError> {
        ni()
    }

    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    pub fn estimate_memory(&self, _s: &Settings) -> Result<MemoryEstimate, CoreError> {
        ni()
    }

    pub fn residual_f64(&self, _field: &Field) -> f64 {
        unimplemented!()
    }

    /// Start a solve session. This consumes the problem.
    pub fn start(self, _s: Settings) -> Result<SolveSession, CoreError> {
        ni()
    }

    /// Solve in one call: a loop over `step`.
    pub fn solve(self, _s: Settings) -> Result<(Field, SolveReport), CoreError> {
        ni()
    }
}
