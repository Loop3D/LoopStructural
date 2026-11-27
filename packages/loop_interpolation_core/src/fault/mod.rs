//! Native fault kernel: restore and apply fault displacement to points.

pub mod gradient;
pub mod grid_field;
pub mod nearest;
#[cfg(feature = "parallel")]
pub mod pool;
pub mod profile;
pub mod region;
pub mod restore;

pub use grid_field::GridField;
pub use profile::{Profile, Profile1};
pub use region::Region;

/// A set of faults that act together on one model.
#[derive(Debug, Clone)]
pub struct FaultChain {
    /// The faults of the chain.
    pub faults: Vec<Fault>,
    /// The settings for the block displacement.
    pub block: BlockSettings,
}

/// One fault with its frame, cutters, regions and displacement profile.
#[derive(Debug, Clone)]
pub struct Fault {
    /// The three frame coordinate fields.
    pub frame: [GridField; 3],
    /// The indices of the faults that cut this fault, in restore order.
    pub cutters: Vec<usize>,
    /// The regions where the fault is active: one list for each frame member.
    /// A member has a value only where all regions of its own list contain the point.
    pub regions: [Vec<Region>; 3],
    /// The displacement profile of the fault.
    pub profile: Profile,
    /// The maximum displacement of the fault.
    pub displacement: f64,
    /// The number of steps in which the displacement is applied.
    pub steps: usize,
}

/// Settings for the displacement of the fault blocks.
#[derive(Debug, Clone)]
pub struct BlockSettings {
    /// The small distance used to decide on which side of the fault a point is.
    pub side_epsilon: f64,
    /// The way in which the profile is continued outside its limits.
    pub extension: Extension,
    /// The value of the frame coordinate inside which the displacement is exact.
    pub exact: f64,
    /// The limit of the frame coordinate for the displacement.
    pub limit: f64,
    /// The step size for the continuation of the profile.
    pub step: f64,
}

/// The way in which a profile is continued outside its limits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Extension {
    /// Continue the profile smoothly.
    Smooth,
    /// Mirror the profile.
    Mirror,
}

/// The side of a fault on which a point is forced to be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// No side is forced (Python 2).
    Free,
    /// The fault is absent (Python 0).
    Absent,
    /// The hanging wall (Python +1).
    Hanging,
    /// The foot wall (Python -1).
    Foot,
}
