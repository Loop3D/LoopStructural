//! Regions in which a fault is active.

use super::restore::force_side;
use super::{FaultChain, GridField, Side};

/// A region of space in which a fault is active.
#[derive(Debug, Clone)]
pub enum Region {
    /// The positive or negative side of a fault.
    Sign {
        /// True for the positive side.
        positive: bool,
        /// The index of the parent fault in the chain.
        fault: usize,
        /// The direction vector of the region, if any.
        vector: Option<[f64; 3]>,
        /// The reference point of the region, if any.
        point: Option<[f64; 3]>,
    },
    /// The side of an unconformity surface.
    Unconformity {
        /// The field of the unconformity.
        field: GridField,
        /// The indices of the faults that restore the point, in restore order.
        faults: Vec<usize>,
        /// The field value of the unconformity surface.
        value: f64,
        /// True if the region has field values at or below `value`.
        sign: bool,
    },
}

impl Region {
    /// Return true if the point `p` is in the region.
    pub fn contains(
        &self,
        chain: &FaultChain,
        p: [f64; 3],
        sides: &[Side],
        restored: &[bool],
    ) -> bool {
        match self {
            Region::Sign {
                positive,
                fault,
                vector,
                point,
            } => {
                // The builder always sets both. Without them the region has no side.
                let (Some(v), Some(c)) = (vector, point) else {
                    return false;
                };
                // A sign region does not take `restored`: the parent is read with no fault restored.
                let none = vec![false; restored.len()];
                let val = force_side(
                    chain.frame_value(*fault, 0, p, sides, &none),
                    sides[*fault],
                    chain.block.side_epsilon,
                );
                let distance = (0..3).map(|i| (c[i] - p[i]) * v[i]).sum::<f64>();
                let distance = force_side(distance, sides[*fault], chain.block.side_epsilon);
                let (val_ok, dist_ok) = if *positive {
                    (val > 0.0, distance > 0.0)
                } else {
                    (val < 0.0, distance < 0.0)
                };
                if val.is_nan() { dist_ok } else { val_ok }
            }
            Region::Unconformity {
                field,
                faults,
                value,
                sign,
            } => {
                let q = chain.restore(faults, p, false, sides, restored);
                let v = field.value(q);
                if *sign { v <= *value } else { v >= *value }
            }
        }
    }
}
