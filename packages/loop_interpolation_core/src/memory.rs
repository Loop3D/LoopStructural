//! Memory tracking and memory estimate.

use crate::error::CoreError;
use crate::report::MemoryEstimate;
use crate::settings::{Precision, Settings};

/// Counts allocated bytes, keeps the peak, and refuses to go over the budget.
#[derive(Debug, Clone)]
pub struct MemTracker {
    current: u64,
    peak: u64,
    #[allow(dead_code)]
    budget: Option<u64>,
}

impl MemTracker {
    pub fn new(budget: Option<u64>) -> Self {
        Self {
            current: 0,
            peak: 0,
            budget,
        }
    }

    /// Allocate a vector of `n` items, each set to `T::default()`.
    ///
    /// The byte count of the request is checked against the budget first.
    /// Nothing is allocated when the request is refused.
    pub fn try_zeroed<T: Copy + Default>(&mut self, n: usize) -> Result<Vec<T>, CoreError> {
        // Compute the request in `u64`. A product that does not fit is out of memory.
        let item_bytes = std::mem::size_of::<T>() as u64;
        let request = (n as u64)
            .checked_mul(item_bytes)
            .ok_or(CoreError::OutOfMemory {
                requested_bytes: u64::MAX,
            })?;

        // Check the budget before any allocation.
        if let Some(budget) = self.budget {
            let estimate = self.current.saturating_add(request);
            if estimate > budget {
                return Err(CoreError::OverBudget {
                    estimate_bytes: estimate,
                    budget_bytes: budget,
                    // The caller fills in the real grid size later.
                    max_nodes_per_axis: 0,
                });
            }
        }

        // Reserve the memory. A failure is returned as an error, not an abort.
        let mut values: Vec<T> = Vec::new();
        values
            .try_reserve_exact(n)
            .map_err(|_| CoreError::OutOfMemory {
                requested_bytes: request,
            })?;
        values.resize(n, T::default());

        self.current = self.current.saturating_add(request);
        self.peak = self.peak.max(self.current);
        Ok(values)
    }

    /// Give back `bytes` to the budget.
    pub fn release(&mut self, bytes: u64) {
        self.current = self.current.saturating_sub(bytes);
    }

    pub fn current(&self) -> u64 {
        self.current
    }

    pub fn peak(&self) -> u64 {
        self.peak
    }
}

/// Counts of a problem. The estimate uses these.
#[derive(Debug, Clone, Default)]
pub struct ProblemCounts {
    /// Number of unknowns (nodes of the grid).
    pub dof: u64,
    pub n_value: u64,
    pub n_gradient: u64,
    pub n_normal: u64,
    pub n_tangent: u64,
    /// Interface points.
    pub n_interface: u64,
    /// Rows of pair bounds `lower <= f(a) - f(b) <= upper`.
    pub n_pair_bound: u64,
    /// Rows of value bounds.
    pub n_value_bound: u64,
    /// True when directional regularisation is on.
    pub uses_directional: bool,
}

// Bytes per unknown (`dof`), with `f64` vectors. A margin of 10 percent
// is added to the sum at the end.
//
// | Part        | Item                                   | Bytes per dof |
// |-------------|----------------------------------------|---------------|
// | vectors     | x, r, p, ap, inv_diag (5 x 8)          | 40            |
// | vectors     | border scratch (1 x 8)                 | 8             |
// | working     | shell scratch (6 slots, small share)   | 8             |
// | working     | shell lookup (u32)                     | 4             |
// | working     | shell row and output indices (u32)     | 3             |
// | working     | copy of the field for the caller (f64) | 8             |
// | working     | ADMM extra vector (none; rhs is rebuilt)| 0             |
// | working     | directional weights (only if used)     | 48            |
//
// Plain CG: 48 + 8 + 4 + 3 + 8 = 71 bytes per dof, 78 with the margin.
// With `f32` the vector part is half (24).
const VECTOR_BYTES_PER_DOF: u64 = 48;
const WORKING_BYTES_PER_DOF: u64 = 8 + 4 + 3 + 8;
const ADMM_BYTES_PER_DOF: u64 = 0;
const DIRECTIONAL_BYTES_PER_DOF: u64 = 48;

// Bytes per constraint record (point, vector or bounds, weight, and the
// row index that the solver keeps).
//
// | Family        | Bytes per row |
// |---------------|---------------|
// | value         | 144           |
// | gradient      | 144           |
// | normal        | 144           |
// | tangent       | 144           |
// | interface     | 288           |
// | pair bound    | 304           |
// | value bound   | 160           |
//
// A stored row of the interpolator holds 8 `usize` indices, 8 `f64`
// coefficients, a value and a weight: 144 bytes. A pair row holds two such
// parts. The counts are stored rows. A normal makes 3 rows.
const VALUE_RECORD: u64 = 144;
const VECTOR_RECORD: u64 = 144;
const INTERFACE_RECORD: u64 = 288;
const PAIR_BOUND_RECORD: u64 = 304;
const VALUE_BOUND_RECORD: u64 = 160;
/// ADMM keeps 3 `f64` values for each bound row (z, dual, work).
const ADMM_BYTES_PER_BOUND_ROW: u64 = 24;

/// Estimate the memory of a solve. It allocates nothing.
pub fn estimate(c: &ProblemCounts, s: &Settings) -> MemoryEstimate {
    let dof = c.dof;

    // Vectors. `f32` halves the vector part.
    let mut vectors = dof.saturating_mul(VECTOR_BYTES_PER_DOF);
    if s.precision == Precision::F32 {
        vectors /= 2;
    }

    // Constraint records.
    let vector_rows = c
        .n_gradient
        .saturating_add(c.n_normal)
        .saturating_add(c.n_tangent);
    let constraints = c
        .n_value
        .saturating_mul(VALUE_RECORD)
        .saturating_add(vector_rows.saturating_mul(VECTOR_RECORD))
        .saturating_add(c.n_interface.saturating_mul(INTERFACE_RECORD))
        .saturating_add(c.n_pair_bound.saturating_mul(PAIR_BOUND_RECORD))
        .saturating_add(c.n_value_bound.saturating_mul(VALUE_BOUND_RECORD));

    // Working memory: tables, scratch, ADMM extra, directional weights.
    let mut per_dof = WORKING_BYTES_PER_DOF;
    let bound_rows = c.n_pair_bound.saturating_add(c.n_value_bound);
    let mut working = 0u64;
    if bound_rows > 0 {
        per_dof += ADMM_BYTES_PER_DOF;
        working = bound_rows.saturating_mul(ADMM_BYTES_PER_BOUND_ROW);
    }
    if c.uses_directional {
        per_dof += DIRECTIONAL_BYTES_PER_DOF;
    }
    working = working.saturating_add(dof.saturating_mul(per_dof));

    // Total with a margin of 10 percent. Use `u128` so that it cannot overflow.
    let sum = vectors as u128 + constraints as u128 + working as u128;
    let total = (sum * 110 / 100).min(u64::MAX as u128) as u64;

    MemoryEstimate {
        vectors,
        constraints,
        working,
        total,
    }
}

/// Largest number of nodes for each axis that fits in `budget`.
/// The result is at least 1. The function allocates nothing.
pub fn max_grid_for_budget(c: &ProblemCounts, s: &Settings, budget: u64) -> usize {
    // 2 000 000 ^ 3 = 8e18, this fits in `u64`.
    const MAX_AXIS: u64 = 2_000_000;
    let fits = |n: u64| {
        let mut counts = c.clone();
        counts.dof = n * n * n;
        estimate(&counts, s).total <= budget
    };
    // Invariant: `lo` fits (or is 1), `hi` is the upper limit.
    let mut lo = 1u64;
    let mut hi = MAX_AXIS;
    while lo < hi {
        let mid = lo + (hi - lo).div_ceil(2);
        if fits(mid) {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    lo as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::CoreError;

    /// Size of one `f64` in bytes.
    const F64_BYTES: u64 = 8;

    // A request for 10 `f64` values counts 80 bytes.
    // The peak becomes 80 too.
    #[test]
    fn allocation_counts_bytes_and_peak() {
        let mut tracker = MemTracker::new(None);
        let v: Vec<f64> = tracker
            .try_zeroed(10)
            .expect("allocation of 10 f64 must succeed");
        assert_eq!(v.len(), 10);
        assert_eq!(tracker.current(), 10 * F64_BYTES);
        assert_eq!(tracker.peak(), 10 * F64_BYTES);
    }

    // Releasing the 80 bytes sets the current count to 0.
    // The peak stays at 80.
    #[test]
    fn release_lowers_current_and_keeps_peak() {
        let mut tracker = MemTracker::new(None);
        let _v: Vec<f64> = tracker
            .try_zeroed(10)
            .expect("allocation of 10 f64 must succeed");
        tracker.release(10 * F64_BYTES);
        assert_eq!(tracker.current(), 0);
        assert_eq!(tracker.peak(), 10 * F64_BYTES);
    }

    // The budget is 79 bytes. A request for 80 bytes is refused.
    // Nothing is counted, so current and peak stay at 0.
    #[test]
    fn over_budget_request_is_refused_before_allocation() {
        let mut tracker = MemTracker::new(Some(79));
        let result = tracker.try_zeroed::<f64>(10);
        assert!(
            matches!(result, Err(CoreError::OverBudget { .. })),
            "expected OverBudget error"
        );
        assert_eq!(tracker.current(), 0);
        assert_eq!(tracker.peak(), 0);
    }

    // A request that is far too large returns an error.
    // It must not abort the process.
    #[test]
    fn huge_request_gives_out_of_memory_not_abort() {
        let mut tracker = MemTracker::new(None);
        let result = tracker.try_zeroed::<f64>(usize::MAX / 16);
        assert!(
            matches!(result, Err(CoreError::OutOfMemory { .. })),
            "expected OutOfMemory error"
        );
    }
}
