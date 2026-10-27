//! Progress and result reports of a solve.

/// Stage of a solve.
#[derive(Debug, Clone)]
pub enum Stage {
    Setup,
    Solve,
    CheckConstraints,
}

/// Progress of a running solve.
#[derive(Debug, Clone)]
pub struct Progress {
    pub stage: Stage,
    pub iteration: usize,
    pub residual: f64,
    pub max_violation: f64,
}

/// Result report of a solve.
#[derive(Debug, Clone)]
pub struct SolveReport {
    pub converged: bool,
    pub iterations: usize,
    pub residual: f64,
    pub max_violation: f64,
    pub peak_bytes: u64,
    pub warnings: Vec<String>,
}

/// Result of one `step` call.
#[derive(Debug, Clone)]
pub enum Step {
    Running(Progress),
    Done(SolveReport),
}

/// Memory estimate in bytes.
#[derive(Debug, Clone)]
pub struct MemoryEstimate {
    pub vectors: u64,
    pub constraints: u64,
    pub working: u64,
    pub total: u64,
}
