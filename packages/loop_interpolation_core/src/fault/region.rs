//! Regions in which a fault is active.

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
        _chain: &FaultChain,
        _p: [f64; 3],
        _sides: &[Side],
        _restored: &[bool],
    ) -> bool {
        unimplemented!("T012")
    }
}
