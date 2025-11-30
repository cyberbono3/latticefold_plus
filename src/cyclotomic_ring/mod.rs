//! Core mathematical infrastructure for LatticeFold+ cyclotomic ring arithmetic
//!
//! Implements the cyclotomic ring `R = Z[X]/(X^d + 1)` with optimized
//! polynomial arithmetic, modular reduction, and SIMD-accelerated operations.
//! The code follows the LatticeFold+ paper for power-of-two dimensions,
//! balanced coefficient representations, and overflow-aware arithmetic.

pub mod constants;
pub mod params;
pub mod balanced_coefficients;
pub mod ring_element;
pub mod ring_config;
pub mod traits;
pub mod error;
pub mod ntt_params;

#[cfg(test)]
mod tests;

pub use balanced_coefficients::BalancedCoefficients;
pub use constants::{MAX_RING_DIMENSION, MEMORY_ALIGNMENT, MIN_RING_DIMENSION, SIMD_WIDTH};
pub use params::RingParams;
pub use ring_element::RingElement;
pub use ring_config::RingConfig;
pub use traits::{NegacyclicMul, RingLike};
pub use error::{CyclotomicRingError, Result};
pub use ntt_params::NTTParams;
