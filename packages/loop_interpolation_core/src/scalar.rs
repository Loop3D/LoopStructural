//! Scalar trait for the vector element type (f32 or f64).
//!
//! Stencil coefficients stay f64. Each sum runs in f64. Only the stored
//! vector values use the type `Scalar`.

use std::ops::{Add, AddAssign, Mul, Sub};

/// Element type of a solver vector. Implemented for `f32` and `f64`.
pub trait Scalar:
    Copy
    + Send
    + Sync
    + Default
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + AddAssign
{
    /// Convert from f64 (rounds for f32).
    fn from_f64(x: f64) -> Self;
    /// Convert to f64 (exact).
    fn to_f64(self) -> f64;
}

impl Scalar for f64 {
    #[inline]
    fn from_f64(x: f64) -> Self {
        x
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self
    }
}

impl Scalar for f32 {
    #[inline]
    fn from_f64(x: f64) -> Self {
        x as f32
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_round_trip_is_exact() {
        let v = 0.123456789012345_f64;
        assert_eq!(f64::from_f64(v).to_f64(), v);
    }

    #[test]
    fn f32_round_trip_is_close() {
        let v = 0.123456789_f64;
        let r = f32::from_f64(v).to_f64();
        assert!((r - v).abs() / v < 1e-7);
        assert_eq!(f32::from_f64(0.5).to_f64(), 0.5);
    }
}
