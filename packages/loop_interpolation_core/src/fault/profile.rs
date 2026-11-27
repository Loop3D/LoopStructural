//! Displacement profiles of a fault.

use super::{BlockSettings, Side};

/// A one-dimensional displacement profile.
#[derive(Debug, Clone)]
pub enum Profile1 {
    /// The value is always one.
    Ones,
    /// The value is always zero.
    Zeros,
    /// A smooth peak function.
    SmoothPeak,
    /// A cubic function with limits.
    Cubic {
        /// The cubic coefficients.
        w: [f64; 4],
        /// The lower and upper limit of the input.
        lim: [f64; 2],
        /// The minimum output value.
        min: f64,
        /// The maximum output value.
        max: f64,
        /// The points (input, value) that define the function.
        value_points: Vec<(f64, f64)>,
    },
    /// One profile for positive input and one for negative input.
    Composite {
        /// The profile for positive input.
        positive: Box<Profile1>,
        /// The profile for negative input.
        negative: Box<Profile1>,
    },
}

impl Profile1 {
    /// Return the value of the profile at `v`.
    pub fn eval(&self, _v: f64) -> f64 {
        unimplemented!("T011")
    }
}

/// The displacement profile of a fault in the three frame coordinates.
#[derive(Debug, Clone)]
pub struct Profile {
    /// The scale of the displacement.
    pub scale: f64,
    /// The profile along frame coordinate 0.
    pub gx: Profile1,
    /// The profile along frame coordinate 1.
    pub gy: Profile1,
    /// The profile along frame coordinate 2.
    pub gz: Profile1,
}

impl Profile {
    /// Return the displacement for the frame coordinates `gx`, `gy` and `gz`.
    pub fn eval(&self, _gx: f64, _gy: f64, _gz: f64) -> f64 {
        unimplemented!("T011")
    }

    /// Return the displacement of a block point, with `forced_gx` as the forced
    /// value of frame coordinate 0 for the given `side`.
    pub fn block_displacement(
        &self,
        _gx: f64,
        _forced_gx: f64,
        _gy: f64,
        _gz: f64,
        _side: Side,
        _block: &BlockSettings,
    ) -> f64 {
        unimplemented!("T011")
    }
}
