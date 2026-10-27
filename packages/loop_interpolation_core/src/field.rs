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

/// Check the length of the points and that each point is inside the grid.
/// Return the points as 3-value arrays.
fn check_points(grid: &StructuredGrid3D, points: &[f64]) -> Result<Vec<[f64; 3]>, CoreError> {
    if !points.len().is_multiple_of(3) {
        return Err(CoreError::InvalidInput {
            constraint: "points".into(),
            reason: format!("length {} is not a multiple of 3", points.len()),
        });
    }
    let pts: Vec<[f64; 3]> = points.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    // The grid code computes (p - origin) / step in f64. NaN is outside.
    for (index, p) in pts.iter().enumerate() {
        if grid.trilinear_coeffs(*p).is_none() {
            return Err(CoreError::PointOutsideGrid { index });
        }
    }
    Ok(pts)
}

impl Field {
    /// Return the values as `f64`.
    pub fn values_f64(&self) -> Vec<f64> {
        match &self.data {
            FieldData::F64(v) => v.clone(),
            FieldData::F32(v) => v.iter().map(|x| *x as f64).collect(),
        }
    }

    /// Evaluate the field. `points` has 3m values. The result has m values.
    pub fn evaluate(&self, points: &[f64]) -> Result<Vec<f64>, CoreError> {
        let pts = check_points(&self.grid, points)?;
        let owned;
        let c: &[f64] = match &self.data {
            FieldData::F64(v) => v,
            FieldData::F32(_) => {
                owned = self.values_f64();
                &owned
            }
        };
        // Same formula as the ported interpolator evaluate method.
        Ok(pts
            .iter()
            .map(|p| match self.grid.trilinear_coeffs(*p) {
                Some((idx, w)) => idx.iter().zip(w.iter()).map(|(i, w)| w * c[*i]).sum(),
                None => f64::NAN,
            })
            .collect())
    }

    /// Evaluate the gradient. `points` has 3m values. The result has 3m values.
    pub fn evaluate_gradient(&self, points: &[f64]) -> Result<Vec<f64>, CoreError> {
        let pts = check_points(&self.grid, points)?;
        let owned;
        let c: &[f64] = match &self.data {
            FieldData::F64(v) => v,
            FieldData::F32(_) => {
                owned = self.values_f64();
                &owned
            }
        };
        let mut out = Vec::with_capacity(points.len());
        for p in &pts {
            match self.grid.trilinear_gradient(*p) {
                Some((idx, grad)) => {
                    for axis in 0..3 {
                        out.push(
                            idx.iter()
                                .zip(grad.iter())
                                .map(|(i, g)| c[*i] * g[axis])
                                .sum(),
                        );
                    }
                }
                None => out.extend([f64::NAN; 3]),
            }
        }
        Ok(out)
    }
}
