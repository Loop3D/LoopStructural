//! Gradient of a faulted field.

use super::{FaultChain, GridField, Side};

/// Return the gradient of `field` at `points` after restore through `order`.
/// The second result is the number of points refined at least once.
#[allow(clippy::too_many_arguments)]
pub fn feature_gradient(
    _chain: &FaultChain,
    _field: &GridField,
    _points: &[[f64; 3]],
    _order: &[usize],
    _sides: &[Side],
    _restored: &[bool],
    _jumping: &[bool],
    _scale: f64,
) -> (Vec<[f64; 3]>, usize) {
    unimplemented!("T024")
}
