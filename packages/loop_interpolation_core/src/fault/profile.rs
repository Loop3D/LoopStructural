//! Displacement profiles of a fault.

use super::{BlockSettings, Extension, Side};

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
    pub fn eval(&self, v: f64) -> f64 {
        match self {
            Profile1::Ones => 1.0,
            Profile1::Zeros => 0.0,
            Profile1::SmoothPeak => {
                if (-1.0..=1.0).contains(&v) {
                    v.powi(4) - 2.0 * v.powi(2) + 1.0
                } else {
                    0.0
                }
            }
            Profile1::Cubic {
                w,
                lim,
                min,
                max,
                value_points,
            } => {
                let limit_value = |x: f64| {
                    value_points
                        .iter()
                        .find(|p| p.0 == x)
                        .map_or_else(|| cubic(w, x), |p| p.1)
                };
                let mut e = cubic(w, v);
                if v >= lim[1] {
                    e = limit_value(lim[1]);
                }
                if v <= lim[0] {
                    e = limit_value(lim[0]);
                }
                if e > *max {
                    e = *max;
                }
                if e < *min {
                    e = *min;
                }
                e
            }
            Profile1::Composite { positive, negative } => {
                if v > 0.0 {
                    positive.eval(v)
                } else if v < 0.0 {
                    negative.eval(v)
                } else {
                    0.0
                }
            }
        }
    }
}

fn cubic(w: &[f64; 4], v: f64) -> f64 {
    w[0] * v.powi(3) + w[1] * v.powi(2) + w[2] * v + w[3]
}

/// Limit a distance from the fault: the identity up to `exact`, constant
/// beyond `limit`, and C2 in between.
fn saturate(distance: f64, block: &BlockSettings) -> f64 {
    let length = block.limit - block.exact;
    let t = ((distance - block.exact) / length).clamp(0.0, 1.0);
    distance.min(block.exact) + length * (t - t.powi(3) + 0.5 * t.powi(4))
}

/// Second-order expansion of `f` at the fault, seen from the side `side`.
fn taylor(f: impl Fn(f64) -> f64, x: f64, side: f64, block: &BlockSettings) -> f64 {
    let h = block.step;
    let x0 = side * block.side_epsilon;
    let k = |i: f64| f(x0 + side * i * h);
    let (f0, f1, f2, f3, f4, f5) = (k(0.0), k(1.0), k(2.0), k(3.0), k(4.0), k(5.0));
    let slope = side * (-25.0 * f0 + 48.0 * f1 - 36.0 * f2 + 16.0 * f3 - 3.0 * f4) / (12.0 * h);
    let curvature =
        (45.0 * f0 - 154.0 * f1 + 214.0 * f2 - 156.0 * f3 + 61.0 * f4 - 10.0 * f5) / (12.0 * h * h);
    f0 + slope * (x - x0) + 0.5 * curvature * (x - x0) * (x - x0)
}

/// Continue the profile of the side `side` to the other side of the fault.
fn continue_profile(profile: &Profile1, x: f64, side: f64, block: &BlockSettings) -> f64 {
    match profile {
        Profile1::Ones | Profile1::Zeros => return profile.eval(x),
        _ => {}
    }
    let x = -side * saturate(x.abs(), block);
    if let Profile1::Cubic { w, .. } = profile {
        // The cubic is the profile at the fault unless a limit, a floor or a ceiling acts there.
        let near = [side * block.side_epsilon, side * block.step];
        if near
            .iter()
            .all(|&n| (profile.eval(n) - cubic(w, n)).abs() <= 1e-12)
        {
            return cubic(w, x);
        }
    }
    taylor(|v| profile.eval(v), x, side, block)
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
    pub fn eval(&self, gx: f64, gy: f64, gz: f64) -> f64 {
        self.scale * self.gx.eval(gx) * self.gy.eval(gy) * self.gz.eval(gz)
    }

    /// Return the displacement of a block point, with `forced_gx` as the forced
    /// value of frame coordinate 0 for the given `side`.
    pub fn block_displacement(
        &self,
        gx: f64,
        forced_gx: f64,
        gy: f64,
        gz: f64,
        side: Side,
        block: &BlockSettings,
    ) -> f64 {
        let s = match side {
            Side::Hanging => 1.0,
            Side::Foot => -1.0,
            _ => return self.eval(forced_gx, gy, gz),
        };
        if block.extension != Extension::Smooth || s * gx >= 0.0 {
            return self.eval(forced_gx, gy, gz);
        }
        let part = match &self.gx {
            Profile1::Composite { positive, negative } => {
                if s > 0.0 {
                    positive
                } else {
                    negative
                }
            }
            other => other,
        };
        self.scale * continue_profile(part, gx, s, block) * self.gy.eval(gy) * self.gz.eval(gz)
    }
}
