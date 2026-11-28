//! Error types of the core crate.

use std::fmt;

/// Kind of a bound constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundKind {
    Value,
    Pair,
}

/// A bound constraint that the field violates.
#[derive(Debug, Clone, PartialEq)]
pub struct ViolatedConstraint {
    pub index: usize,
    pub kind: BoundKind,
    pub violation: f64,
}

/// Error of the core crate.
#[derive(Debug, Clone, PartialEq)]
pub enum CoreError {
    /// A constraint is not valid. `constraint` is a text such as "value[12]".
    InvalidInput { constraint: String, reason: String },
    /// An axis has fewer than 3 nodes.
    GridTooSmall { axis: usize, nodes: usize },
    /// The grid has more nodes than `max_nodes`.
    GridTooLarge { nodes: u64, max_nodes: u64 },
    /// The problem has no constraints.
    NoConstraints,
    /// The memory estimate is more than the budget.
    OverBudget {
        estimate_bytes: u64,
        budget_bytes: u64,
        max_nodes_per_axis: usize,
    },
    /// An allocation failed.
    OutOfMemory { requested_bytes: u64 },
    /// The bounds cannot all be true. At most 20 entries, largest violation first.
    Infeasible { worst: Vec<ViolatedConstraint> },
    /// A point is outside the grid.
    PointOutsideGrid { index: usize },
    /// A setting is not supported.
    Unsupported { setting: String },
    /// A solve is already running (used by the wasm layer).
    SolveRunning,
}

impl fmt::Display for CoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CoreError::InvalidInput { constraint, reason } => {
                write!(f, "invalid input in {constraint}: {reason}")
            }
            CoreError::GridTooSmall { axis, nodes } => {
                write!(
                    f,
                    "grid axis {axis} has {nodes} nodes; at least 3 are needed"
                )
            }
            CoreError::GridTooLarge { nodes, max_nodes } => {
                write!(f, "grid has {nodes} nodes; the maximum is {max_nodes}")
            }
            CoreError::NoConstraints => write!(f, "the problem has no constraints"),
            CoreError::OverBudget {
                estimate_bytes,
                budget_bytes,
                max_nodes_per_axis,
            } => write!(
                f,
                "estimated memory {estimate_bytes} bytes is more than the budget {budget_bytes} bytes; the largest grid for the budget has {max_nodes_per_axis} nodes per axis"
            ),
            CoreError::OutOfMemory { requested_bytes } => {
                write!(
                    f,
                    "out of memory: could not allocate {requested_bytes} bytes"
                )
            }
            CoreError::Infeasible { worst } => {
                write!(
                    f,
                    "the bounds are infeasible; {} worst violations listed",
                    worst.len()
                )
            }
            CoreError::PointOutsideGrid { index } => {
                write!(f, "point {index} is outside the grid")
            }
            CoreError::Unsupported { setting } => write!(f, "unsupported setting: {setting}"),
            CoreError::SolveRunning => write!(f, "a solve is already running"),
        }
    }
}

impl std::error::Error for CoreError {}
