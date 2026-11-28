//! Test helper: a counting allocator.
//! It wraps the system allocator and tracks the bytes now in use (CURRENT)
//! and the highest value of CURRENT since the last reset (PEAK).
#![allow(dead_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct CountingAllocator;

static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

fn add(bytes: usize) {
    let now = CURRENT.fetch_add(bytes, Ordering::SeqCst) + bytes;
    PEAK.fetch_max(now, Ordering::SeqCst);
}

fn sub(bytes: usize) {
    CURRENT.fetch_sub(bytes, Ordering::SeqCst);
}

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(layout) };
        if !p.is_null() {
            add(layout.size());
        }
        p
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { System.alloc_zeroed(layout) };
        if !p.is_null() {
            add(layout.size());
        }
        p
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
        sub(layout.size());
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let p = unsafe { System.realloc(ptr, layout, new_size) };
        if !p.is_null() {
            // On success the old block is gone and the new block exists.
            if new_size >= layout.size() {
                add(new_size - layout.size());
            } else {
                sub(layout.size() - new_size);
            }
        }
        p
    }
}

#[global_allocator]
static ALLOC: CountingAllocator = CountingAllocator;

/// Set PEAK to the bytes now in use.
pub fn reset_peak() {
    PEAK.store(CURRENT.load(Ordering::SeqCst), Ordering::SeqCst);
}

/// Bytes now in use.
pub fn current() -> usize {
    CURRENT.load(Ordering::SeqCst)
}

/// Highest bytes in use since the last reset.
pub fn peak() -> usize {
    PEAK.load(Ordering::SeqCst)
}
