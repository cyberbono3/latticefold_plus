//! Core mathematical infrastructure for LatticeFold+ cyclotomic ring arithmetic
//!
//! Implements the cyclotomic ring `R = Z[X]/(X^d + 1)` with optimized
//! polynomial arithmetic, modular reduction, and SIMD-accelerated operations.
//! The code follows the LatticeFold+ paper for power-of-two dimensions,
//! balanced coefficient representations, and overflow-aware arithmetic.

pub mod constants;
pub mod error;

pub mod params;
pub mod ring;
pub mod traits;

pub use constants::{
    MAX_RING_DIMENSION, MEMORY_ALIGNMENT, MIN_RING_DIMENSION, SIMD_WIDTH,
};
pub use error::{CyclotomicRingError, Result};
pub use params::RingParams;
pub use ring::CyclotomicRing;
pub use traits::Ring;
