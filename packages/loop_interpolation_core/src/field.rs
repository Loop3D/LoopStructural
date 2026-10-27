//! The solved field.

use crate::error::CoreError;
use crate::grid::StructuredGrid3D;

/// Values of the field, in the solver precision.
#[derive(Debug, Clone)]
pub enum FieldData {
    F64(Vec<f64>),
    F32(Vec<f32>),
}

/// A solved field and its grid.
#[derive(Debug, Clone)]
pub struct Field {
    pub grid: StructuredGrid3D,
    pub data: FieldData,
}

fn not_implemented() -> CoreError {
    CoreError::Unsupported {
        setting: "not implemented".into(),
    }
}

impl Field {
    /// Return the values as `f64`.
    pub fn values_f64(&self) -> Vec<f64> {
        unimplemented!()
    }

    /// Evaluate the field. `points` has 3m values. The result has m values.
    pub fn evaluate(&self, _points: &[f64]) -> Result<Vec<f64>, CoreError> {
        Err(not_implemented())
    }

    /// Evaluate the gradient. `points` has 3m values. The result has 3m values.
    pub fn evaluate_gradient(&self, _points: &[f64]) -> Result<Vec<f64>, CoreError> {
        Err(not_implemented())
    }
}
