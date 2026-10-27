//! Solver settings.

/// Floating-point type of the solver vectors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Precision {
    F64,
    F32,
}

/// Settings of a solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Settings {
    pub precision: Precision,
    /// Relative residual tolerance.
    pub tol: f64,
    /// Total inner CG iterations.
    pub max_iter: usize,
    pub budget_bytes: Option<u64>,
    /// ADMM penalty.
    pub rho: f64,
    pub inner_max_iter: usize,
    pub inner_tol: f64,
    pub admm_abs_tol: f64,
    pub admm_rel_tol: f64,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            precision: Precision::F64,
            tol: 1e-8,
            max_iter: 5000,
            budget_bytes: None,
            rho: 1.0,
            inner_max_iter: 3000,
            inner_tol: 1e-6,
            admm_abs_tol: 1e-4,
            admm_rel_tol: 1e-3,
        }
    }
}
