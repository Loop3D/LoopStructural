//! Core crate for the low-memory interpolator.
//!
//! The crate has a matrix-free finite-difference interpolator. It does not
//! build a sparse matrix. It applies the constraints as linear operators and
//! solves with a conjugate gradient method.

pub mod cg;
pub mod grid;
pub mod interpolator;
pub(crate) mod par;
pub mod stencil;

/// Rayon chunk-splitting floor shared by every parallel pass in the CG hot
/// loop (`cg.rs`, `interpolator.rs`, `stencil.rs`). Rayon's default
/// splitting for a plain `.par_iter()` recursively halves a slice down to
/// very small chunks (adaptively work-stolen), which -- confirmed by
/// benchmarking the real-world `examples/claudius_fdi_comparison.py`
/// dataset (dof in the tens of thousands, well below the `dof=1_000_000`
/// grid `interpolator.rs`'s `profile_apply_1e6` tests use to validate this
/// crate's performance) -- creates *more* tiny tasks than the ~20 hardware
/// threads can usefully absorb, so thread-pool dispatch/work-stealing
/// overhead (paid on a dozen-plus separate `.par_iter()` calls per CG
/// iteration) dominates over the arithmetic. Forcing complete
/// sequential execution below some dof threshold turned out to be the wrong
/// fix, though: at dof ~45,000 it was *worse* than plain rayon defaults
/// (measured ~7.8s vs. ~3.3s), since that grid has enough real per-node work
/// to profit from spreading across cores -- it just needs coarser, not
/// finer, chunks. `.with_min_len(PAR_MIN_LEN)` caps how many chunks rayon is
/// allowed to split a call into (roughly `len / PAR_MIN_LEN` at most),
/// trading fewer, larger chunks (less dispatch overhead, still parallel)
/// for the very fine adaptive splitting rayon defaults to. 200 was chosen by
/// sweeping candidates against the claudius benchmark at both dof~6,000 and
/// dof~45,000 and picking the value that won at both scales (~3x faster at
/// the small end, ~1.8x at the large end, vs. no `with_min_len` at all) --
/// re-sweep (`examples/claudius_fdi_comparison.py`, plus
/// `profile_apply_1e6`/`profile_apply_1e6_unfused` for the large-grid end)
/// if the target hardware's core count changes significantly.
pub(crate) const PAR_MIN_LEN: usize = 200;

pub use cg::CgResult;
pub use grid::StructuredGrid3D;
pub use interpolator::{AdmmResult, FiniteDifferenceInterpolator, RegularisationWeights};
