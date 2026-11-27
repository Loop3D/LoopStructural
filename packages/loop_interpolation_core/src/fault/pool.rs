//! Thread pools for the fault kernel.

use std::sync::Arc;

/// Return the thread pool for `threads` threads. Zero means all cores.
pub fn pool(_threads: usize) -> Arc<rayon::ThreadPool> {
    unimplemented!("T012")
}
