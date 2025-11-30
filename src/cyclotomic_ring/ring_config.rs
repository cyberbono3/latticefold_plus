use super::error::{CyclotomicRingError, Result};

use super::constants::{MAX_RING_DIMENSION, MIN_RING_DIMENSION};

/// Centralized validation and balanced-reduction utilities for ring arithmetic.
pub struct RingConfig;

impl RingConfig {
    /// Ensures modulus is positive and returns half modulus used for bounds.
    pub fn validate_modulus(modulus: i64) -> Result<i64> {
        if modulus <= 0 {
            Err(CyclotomicRingError::InvalidModulus { modulus })
        } else {
            Ok(modulus / 2)
        }
    }

    /// Ensures dimension is a supported power of two.
    pub fn validate_dimension(dimension: usize) -> Result<()> {
        if !dimension.is_power_of_two() {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: dimension.next_power_of_two(),
                got: dimension,
            });
        }

        if dimension < MIN_RING_DIMENSION || dimension > MAX_RING_DIMENSION {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: MIN_RING_DIMENSION,
                got: dimension,
            });
        }

        Ok(())
    }

    /// Reduces a single coefficient into balanced representation for modulus q.
    pub fn reduce_to_balanced(coeff: i64, modulus: i64) -> i64 {
        let half_modulus = modulus / 2;
        let mut reduced = coeff % modulus;
        if reduced > half_modulus {
            reduced - modulus
        } else if reduced < -half_modulus {
            reduced + modulus
        } else {
            reduced
        }
    }

    /// Reduces a slice of coefficients into balanced representation for modulus q.
    pub fn reduce_slice_to_balanced(coeffs: &mut [i64], modulus: i64) {
        for coeff in coeffs {
            *coeff = Self::reduce_to_balanced(*coeff, modulus);
        }
    }
}
