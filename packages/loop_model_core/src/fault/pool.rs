//! Thread pools for the fault kernel.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

static POOLS: OnceLock<Mutex<HashMap<usize, Arc<rayon::ThreadPool>>>> = OnceLock::new();

/// Return the thread pool for `threads` threads. Zero means all cores.
/// The function keeps one pool for each thread count.
pub fn pool(threads: usize) -> Arc<rayon::ThreadPool> {
    let mut pools = POOLS
        .get_or_init(Default::default)
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    pools
        .entry(threads)
        .or_insert_with(|| {
            Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("cannot build the thread pool"),
            )
        })
        .clone()
}
