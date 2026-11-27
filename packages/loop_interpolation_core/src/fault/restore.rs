//! Apply and restore fault displacement.

use super::{FaultChain, Side};

/// Put the frame value `v` on the side `side` of the fault. NaN stays NaN.
/// A side other than the hanging wall or the foot wall gives `v` unchanged.
pub fn force_side(v: f64, side: Side, epsilon: f64) -> f64 {
    let s = match side {
        Side::Hanging => 1.0,
        Side::Foot => -1.0,
        _ => return v,
    };
    if v.is_nan() {
        v
    } else {
        s * v.abs().max(epsilon)
    }
}

impl FaultChain {
    /// Move the point `p` through one fault. Use `reverse` to undo the displacement.
    pub fn apply(
        &self,
        fault: usize,
        p: [f64; 3],
        reverse: bool,
        sides: &[Side],
        restored: &[bool],
    ) -> [f64; 3] {
        if sides[fault] == Side::Absent {
            return p;
        }
        let f = &self.faults[fault];
        let sign = if reverse { -1.0 } else { 1.0 };
        let [mut gx, mut gy, mut gz] = self.frame_values(fault, p, sides, restored);
        let fraction = self.displacement_fraction(fault, [gx, gy, gz], sides);
        // The mask uses the fraction before the displacement scale, as the Python code does.
        if fraction.abs().partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return p;
        }
        let mut p = p;
        for step in 0..f.steps {
            if step > 0 {
                [gx, gy, gz] = self.frame_values(fault, p, sides, restored);
            }
            let q = self.frame_point(fault, p, sides, restored);
            let mut g = f.frame[1].gradient(q);
            // A cut frame has no single gradient on a cutting fault (#425): use the plain grid gradient.
            if g.iter().any(|v| v.is_nan())
                && !p.iter().any(|v| v.is_nan())
                && !f.cutters.is_empty()
            {
                g = f.frame[1].gradient(p);
            }
            let ok = !(gx.is_nan() || gy.is_nan() || gz.is_nan());
            let d = self.displacement_fraction(fault, [gx, gy, gz], sides) * f.displacement * sign;
            let mag = if ok {
                g.iter().map(|v| v * v).sum::<f64>().sqrt()
            } else {
                0.0
            };
            if mag > 0.0 {
                g = g.map(|v| v / mag);
            }
            let scale = (1.0 / f.steps as f64) * d;
            for i in 0..3 {
                p[i] += g[i] * scale;
            }
        }
        p
    }

    /// Move the point `p` through the faults in `order`. A fault that is already
    /// restored is skipped, and each fault is marked done after it moves the point.
    pub fn restore(
        &self,
        order: &[usize],
        p: [f64; 3],
        reverse: bool,
        sides: &[Side],
        restored: &[bool],
    ) -> [f64; 3] {
        let mut done = restored.to_vec();
        let mut p = p;
        for &f in order {
            if done[f] {
                continue;
            }
            p = self.apply(f, p, reverse, sides, &done);
            done[f] = true;
        }
        p
    }

    /// Return one frame coordinate of a fault at `p`, with the region mask applied.
    pub fn frame_value(
        &self,
        fault: usize,
        coordinate: usize,
        p: [f64; 3],
        sides: &[Side],
        restored: &[bool],
    ) -> f64 {
        let f = &self.faults[fault];
        if !f.regions[coordinate]
            .iter()
            .all(|r| r.contains(self, p, sides, restored))
        {
            return f64::NAN;
        }
        f.frame[coordinate].value(self.frame_point(fault, p, sides, restored))
    }

    /// Move all points through the faults in `order`. Use `threads` threads
    /// (0 is all cores) with the feature `parallel`; otherwise `threads` has no effect.
    pub fn restore_batch(
        &self,
        order: &[usize],
        points: &[[f64; 3]],
        reverse: bool,
        sides: &[Side],
        restored: &[bool],
        threads: usize,
    ) -> Vec<[f64; 3]> {
        self.batch(points, threads, |p| {
            self.restore(order, p, reverse, sides, restored)
        })
    }

    /// Return one frame coordinate of a fault at all points.
    pub fn frame_value_batch(
        &self,
        fault: usize,
        coordinate: usize,
        points: &[[f64; 3]],
        sides: &[Side],
        restored: &[bool],
        threads: usize,
    ) -> Vec<f64> {
        self.batch(points, threads, |p| {
            self.frame_value(fault, coordinate, p, sides, restored)
        })
    }

    #[cfg(feature = "parallel")]
    fn batch<T: Send>(
        &self,
        points: &[[f64; 3]],
        threads: usize,
        f: impl Fn([f64; 3]) -> T + Sync + Send,
    ) -> Vec<T> {
        use rayon::prelude::*;
        super::pool::pool(threads).install(|| points.par_iter().map(|&p| f(p)).collect())
    }

    #[cfg(not(feature = "parallel"))]
    fn batch<T>(&self, points: &[[f64; 3]], _threads: usize, f: impl Fn([f64; 3]) -> T) -> Vec<T> {
        points.iter().map(|&p| f(p)).collect()
    }

    /// Return the point at which the frame of `fault` is read: `p` restored
    /// through the cutters that are not yet restored.
    fn frame_point(
        &self,
        fault: usize,
        p: [f64; 3],
        sides: &[Side],
        restored: &[bool],
    ) -> [f64; 3] {
        let cutters = &self.faults[fault].cutters;
        if cutters.iter().all(|&c| restored[c]) {
            p
        } else {
            self.restore(cutters, p, false, sides, restored)
        }
    }

    /// Return the three frame coordinates of `fault` at `p`.
    fn frame_values(
        &self,
        fault: usize,
        p: [f64; 3],
        sides: &[Side],
        restored: &[bool],
    ) -> [f64; 3] {
        std::array::from_fn(|i| self.frame_value(fault, i, p, sides, restored))
    }

    /// Return the displacement at the frame coordinates as a fraction of the displacement.
    /// Return zero if a coordinate is NaN.
    fn displacement_fraction(&self, fault: usize, g: [f64; 3], sides: &[Side]) -> f64 {
        let side = sides[fault];
        let forced = force_side(g[0], side, self.block.side_epsilon);
        if forced.is_nan() || g[1].is_nan() || g[2].is_nan() {
            return 0.0;
        }
        self.faults[fault]
            .profile
            .block_displacement(g[0], forced, g[1], g[2], side, &self.block)
    }
}
