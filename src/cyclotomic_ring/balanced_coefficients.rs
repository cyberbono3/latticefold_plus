use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use zeroize::{Zeroize, ZeroizeOnDrop};

use super::error::{CyclotomicRingError, Result};

use super::constants::{MEMORY_ALIGNMENT, SIMD_WIDTH};
use super::ring_config::RingConfig;

/// Balanced coefficients held in the symmetric set Zq = {-⌊q/2⌋, …, ⌊q/2⌋}.
#[derive(Clone, PartialEq, Eq, Zeroize, ZeroizeOnDrop)]
pub struct BalancedCoefficients {
    /// Coefficient vector in balanced representation: coeffs[i] ∈ [-⌊q/2⌋, ⌊q/2⌋]
    pub(crate) coeffs: Vec<i64>,
    
    /// Modulus q for the coefficient ring Zq
    pub(crate) modulus: i64,
    
    /// Precomputed value ⌊q/2⌋ for efficient bounds checking
    pub(crate) half_modulus: i64,
}

impl BalancedCoefficients {
    /// Ensures coefficients are inside the balanced interval [-half_modulus, half_modulus].
    fn validate_balanced_range(coeffs: &[i64], half_modulus: i64) -> Result<()> {
        for (i, &coeff) in coeffs.iter().enumerate() {
            if coeff < -half_modulus || coeff > half_modulus {
                return Err(CyclotomicRingError::CoefficientOutOfRange {
                    coefficient: coeff,
                    min_bound: -half_modulus,
                    max_bound: half_modulus,
                    position: i,
                });
            }
        }
        Ok(())
    }

    /// Build balanced coefficients from a user-supplied vector and modulus.
    pub fn new(coeffs: Vec<i64>, modulus: i64) -> Result<Self> {
        let half_modulus = RingConfig::validate_modulus(modulus)?;
        Self::validate_balanced_range(&coeffs, half_modulus)?;

        Ok(Self {
            coeffs,
            modulus,
            half_modulus,
        })
    }
    
    /// Allocate zeroed, balanced coefficients for a given dimension/modulus.
    pub fn with_dimension(dimension: usize, modulus: i64) -> Result<Self> {
        // Validate inputs
        RingConfig::validate_dimension(dimension)?;
        let half_modulus = RingConfig::validate_modulus(modulus)?;
        
        // Calculate padded dimension for SIMD alignment
        // Padding ensures vectorized operations don't access invalid memory
        let padded_dimension = ((dimension + SIMD_WIDTH - 1) / SIMD_WIDTH) * SIMD_WIDTH;
        
        // Allocate coefficient vector with zero initialization
        // Zero is always in the balanced representation range
        let mut coeffs = vec![0i64; padded_dimension];
        
        // Ensure memory alignment for SIMD operations
        // This is critical for performance on modern CPUs with vector units
        if coeffs.as_ptr() as usize % MEMORY_ALIGNMENT != 0 {
            // Reallocate with proper alignment if needed
            // This should be rare with modern allocators but ensures correctness
            let mut aligned_coeffs = Vec::with_capacity(padded_dimension + MEMORY_ALIGNMENT / 8);
            aligned_coeffs.resize(padded_dimension, 0i64);
            coeffs = aligned_coeffs;
        }
        
        // Truncate to actual dimension (remove padding from logical view)
        // Physical memory remains padded but logical operations use correct size
        coeffs.truncate(dimension);
        
        Ok(Self {
            coeffs,
            modulus,
            half_modulus,
        })
    }
    
    /// Convert a standard-representation slice into balanced coefficients.
    pub fn from_standard(standard_coeffs: &[i64], modulus: i64) -> Result<Self> {
        // Validate input parameters
        // Empty coefficient vector is not meaningful for polynomial operations
        if standard_coeffs.is_empty() {
            return Err(CyclotomicRingError::InvalidParameters(
                "Coefficient vector cannot be empty".to_string(),
            ));
        }
        
        let half_modulus = RingConfig::validate_modulus(modulus)?;
        
        // Validate all coefficients are in standard range [0, q-1]
        // This prevents undefined behavior in conversion
        for (i, &coeff) in standard_coeffs.iter().enumerate() {
            if coeff < 0 || coeff >= modulus {
                return Err(CyclotomicRingError::CoefficientOutOfRange {
                    coefficient: coeff,
                    min_bound: 0,
                    max_bound: modulus - 1,
                    position: i,
                });
            }
        }
        
        // Create new instance with proper dimension and modulus
        let mut result = Self::with_dimension(standard_coeffs.len(), modulus)?;
        for (i, &coeff) in standard_coeffs.iter().enumerate() {
            result.coeffs[i] = RingConfig::reduce_to_balanced(coeff, modulus);
        }
        
        Ok(result)
    }
    
    /// Converts balanced coefficients to standard representation
    pub fn to_standard(&self) -> Vec<i64> {
        // Allocate result vector with same dimension
        let mut standard_coeffs = vec![0i64; self.coeffs.len()];

        for (i, &coeff) in self.coeffs.iter().enumerate() {
            standard_coeffs[i] = if coeff < 0 {
                coeff + self.modulus
            } else {
                coeff
            };
        }
        
        standard_coeffs
    }
    
    /// Returns the coefficient vector as a slice
    pub fn coefficients(&self) -> &[i64] {
        &self.coeffs
    }
    
    /// Validates that all coefficients are within balanced representation bounds
    pub fn validate_bounds(&self) -> Result<()> {
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            if coeff < -self.half_modulus || coeff > self.half_modulus {
                return Err(CyclotomicRingError::CoefficientOutOfRange {
                    coefficient: coeff,
                    min_bound: -self.half_modulus,
                    max_bound: self.half_modulus,
                    position: i,
                });
            }
        }
        Ok(())
    }
    
    /// Returns the modulus used for this coefficient representation
    pub fn modulus(&self) -> i64 {
        self.modulus
    }
    
    /// Returns the dimension (number of coefficients)
    pub fn dimension(&self) -> usize {
        self.coeffs.len()
    }
    
    /// Computes the ℓ∞-norm of the coefficient vector
    pub fn infinity_norm(&self) -> i64 {
        // Handle empty coefficient vector
        if self.coeffs.is_empty() {
            return 0;
        }
        
        // Initialize maximum with first coefficient's absolute value
        let mut max_abs = self.coeffs[0].abs();
        
        for &coeff in &self.coeffs {
            max_abs = max_abs.max(coeff.abs());
        }
        
        max_abs
    }
}

impl Debug for BalancedCoefficients {
    /// Custom debug formatting for balanced coefficients
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_struct("BalancedCoefficients")
            .field("dimension", &self.coeffs.len())
            .field("modulus", &self.modulus)
            .field("coefficients", &self.coeffs)
            .field("infinity_norm", &self.infinity_norm())
            .finish()
    }
}

impl Display for BalancedCoefficients {
    /// User-friendly display formatting for balanced coefficients
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "BalancedCoefficients(dim={}, q={}, ||·||_∞={})", 
               self.coeffs.len(), self.modulus, self.infinity_norm())
    }
}
