//! Memory tracking and memory estimate.

use crate::error::CoreError;
use crate::report::MemoryEstimate;
use crate::settings::Settings;

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
#[derive(Debug, Clone)]
pub struct ProblemCounts {
    #[allow(dead_code)]
    pub(crate) dof: u64,
    #[allow(dead_code)]
    pub(crate) uses_directional: bool,
}

/// Estimate the memory of a solve. It allocates nothing.
pub fn estimate(_c: &ProblemCounts, _s: &Settings) -> MemoryEstimate {
    unimplemented!()
}

/// Largest number of nodes for each axis that fits in `budget`.
pub fn max_grid_for_budget(_c: &ProblemCounts, _s: &Settings, _budget: u64) -> usize {
    unimplemented!()
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
