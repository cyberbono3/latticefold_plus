use super::error::{CyclotomicRingError, Result};

use super::constants::{MAX_RING_DIMENSION, MIN_RING_DIMENSION};

/// Parameters for cyclotomic ring operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingParams {
    /// Ring dimension d (must be power of 2, 32 ≤ d ≤ 16384)
    pub dimension: usize,

    /// Modulus q for coefficient ring Zq
    pub modulus: i64,
}

impl RingParams {
    /// Creates new ring parameters with validation.
    pub fn new(dimension: usize, modulus: i64) -> Result<Self> {
        // Validate dimension is power of 2 within supported range
        if !dimension.is_power_of_two() {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: dimension.next_power_of_two(),
                got: dimension,
            });
        }

        if !(MIN_RING_DIMENSION..=MAX_RING_DIMENSION).contains(&dimension) {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: MIN_RING_DIMENSION,
                got: dimension,
            });
        }

        // Validate modulus is positive
        if modulus <= 0 {
            return Err(CyclotomicRingError::InvalidModulus { modulus });
        }

        Ok(Self { dimension, modulus })
    }

    /// Returns standard parameters for testing and development.
    pub fn standard() -> Self {
        Self {
            dimension: 64,
            modulus: 2_147_483_647, // Large prime
        }
    }
}
