//! Apply and restore fault displacement.

use super::{FaultChain, Side};

impl FaultChain {
    /// Move the point `p` through one fault. Use `reverse` to undo the displacement.
    pub fn apply(
        &self,
        _fault: usize,
        _p: [f64; 3],
        _reverse: bool,
        _sides: &[Side],
        _restored: &[bool],
    ) -> [f64; 3] {
        unimplemented!("T012")
    }

    /// Move the point `p` through the faults in `order`.
    pub fn restore(
        &self,
        _order: &[usize],
        _p: [f64; 3],
        _reverse: bool,
        _sides: &[Side],
        _restored: &[bool],
    ) -> [f64; 3] {
        unimplemented!("T012")
    }

    /// Return one frame coordinate of a fault at `p`, with the region mask applied.
    pub fn frame_value(
        &self,
        _fault: usize,
        _coordinate: usize,
        _p: [f64; 3],
        _sides: &[Side],
        _restored: &[bool],
    ) -> f64 {
        unimplemented!("T012")
    }
}
