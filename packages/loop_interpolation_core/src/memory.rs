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

    /// Allocate a zeroed vector of `n` items.
    pub fn try_zeroed<T: Copy + Default>(&mut self, _n: usize) -> Result<Vec<T>, CoreError> {
        Err(CoreError::Unsupported {
            setting: "not implemented".into(),
        })
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
