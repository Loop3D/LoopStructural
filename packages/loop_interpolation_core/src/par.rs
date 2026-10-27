//! Optional parallel layer.
//!
//! With feature `parallel`, this module gives the `rayon` prelude.
//! Without the feature (for example on `wasm32-unknown-unknown`, which has
//! no threads), it gives traits with the same method names on std
//! iterators and slices. Each method runs on one thread and gives the same
//! result as the `rayon` method.
//!
//! The shim has only the methods that this crate uses: `par_iter`,
//! `par_iter_mut`, `par_chunks_mut`, `into_par_iter` and `with_min_len`.
//! All other iterator methods (`zip`, `enumerate`, `map`, `sum`,
//! `collect`, `for_each`) come from `std::iter::Iterator`.

#[cfg(feature = "parallel")]
pub(crate) use rayon::prelude::*;

#[cfg(not(feature = "parallel"))]
pub(crate) use shim::*;

#[cfg(not(feature = "parallel"))]
mod shim {
    /// Slice methods with the `rayon` names.
    pub(crate) trait ParSliceShim<T> {
        fn par_iter(&self) -> std::slice::Iter<'_, T>;
        fn par_iter_mut(&mut self) -> std::slice::IterMut<'_, T>;
        fn par_chunks_mut(&mut self, chunk_size: usize) -> std::slice::ChunksMut<'_, T>;
    }

    impl<T> ParSliceShim<T> for [T] {
        #[inline]
        fn par_iter(&self) -> std::slice::Iter<'_, T> {
            self.iter()
        }

        #[inline]
        fn par_iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
            self.iter_mut()
        }

        #[inline]
        fn par_chunks_mut(&mut self, chunk_size: usize) -> std::slice::ChunksMut<'_, T> {
            self.chunks_mut(chunk_size)
        }
    }

    /// `into_par_iter` for any value that has `into_iter` (for example a range).
    pub(crate) trait IntoParIterShim: IntoIterator + Sized {
        #[inline]
        fn into_par_iter(self) -> Self::IntoIter {
            self.into_iter()
        }
    }

    impl<I: IntoIterator> IntoParIterShim for I {}

    /// `with_min_len` does nothing with one thread. It returns `self`.
    pub(crate) trait WithMinLenShim: Iterator + Sized {
        #[inline]
        fn with_min_len(self, _min: usize) -> Self {
            self
        }
    }

    impl<I: Iterator> WithMinLenShim for I {}
}
