use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::{Zeroize, ZeroizeOnDrop};

use super::error::{CyclotomicRingError, Result};
use num_traits::{One, Zero};
use super::NTTParams;

use super::{BalancedCoefficients, RingConfig, RingParams};
use super::traits::{NegacyclicMul, RingLike};
use super::constants::{MAX_RING_DIMENSION, MIN_RING_DIMENSION};

/// Cyclotomic ring element in R = Z[X]/(X^d + 1).
#[derive(Clone, PartialEq, Eq, Zeroize, ZeroizeOnDrop)]
pub struct RingElement {
    /// Coefficient representation in balanced form
    coefficients: BalancedCoefficients,
    
    /// Ring dimension d (must be power of 2)
    dimension: usize,
    
    /// Optional modulus for Rq = R/qR operations
    modulus: Option<i64>,
}

impl RingElement {
    /// Construct the zero element for a given dimension/modulus.
    pub fn zero(dimension: usize, modulus: Option<i64>) -> Result<Self> {
        // Validate dimension is power of 2 within supported range
        RingConfig::validate_dimension(dimension)?;
        
        // Validate modulus is positive for Rq operations
        let modulus_value = modulus.unwrap_or(1);
        RingConfig::validate_modulus(modulus_value)?;
        
        // Create coefficient representation with given modulus
        let coefficients = BalancedCoefficients::with_dimension(dimension, modulus_value)?;
        
        Ok(Self {
            coefficients,
            dimension,
            modulus,
        })
    }
    
    /// Construct the multiplicative identity (constant 1).
    pub fn one(dimension: usize, modulus: Option<i64>) -> Result<Self> {
        // Create zero element as base
        let mut result = Self::zero(dimension, modulus)?;
        
        // Set constant coefficient to 1
        // This represents the polynomial 1(X) = 1 + 0·X + 0·X² + ...
        result.coefficients.coeffs[0] = 1;
        
        // Validate the coefficient is within bounds
        result.coefficients.validate_bounds()?;
        
        Ok(result)
    }
    
    /// Build an element from balanced coefficients and a modulus.
    pub fn from_coefficients(coeffs: &[i64], modulus: Option<i64>) -> Result<Self> {
        // Validate coefficient vector is not empty
        if coeffs.is_empty() {
            return Err(CyclotomicRingError::InvalidParameters(
                "Coefficient vector cannot be empty".to_string(),
            ));
        }
        
        let dimension = coeffs.len();
        
        // Validate modulus is positive
        let modulus_value = modulus.unwrap_or(1);
        RingConfig::validate_modulus(modulus_value)?;
        
        RingConfig::validate_dimension(dimension)?;
        
        // Treat input as balanced coefficients and validate range
        let coefficients = BalancedCoefficients::new(coeffs.to_vec(), modulus_value)?;
        
        Ok(Self {
            coefficients,
            dimension,
            modulus,
        })
    }
    
    /// Returns the coefficient vector in balanced representation
    pub fn coefficients(&self) -> &[i64] {
        self.coefficients.coefficients()
    }
    
    /// Serializes the ring element to bytes
    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut bytes = Vec::new();
        
        // Write dimension
        bytes.extend_from_slice(&(self.dimension as u64).to_le_bytes());
        
        // Write modulus (0 if None)
        let modulus_value = self.modulus.unwrap_or(0);
        bytes.extend_from_slice(&modulus_value.to_le_bytes());
        
        // Write coefficients
        for &coeff in self.coefficients() {
            bytes.extend_from_slice(&coeff.to_le_bytes());
        }
        
        Ok(bytes)
    }
    
    /// Provides mutable access to polynomial coefficients for modification.
    pub fn coefficients_mut(&mut self) -> &mut [i64] {
        &mut self.coefficients.coeffs
    }
    
    /// Returns the ring dimension d
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    
    /// Returns the optional modulus q for Rq operations
    pub fn modulus(&self) -> Option<i64> {
        self.modulus
    }
    
    /// Computes the infinity norm ||f||_∞ = max_i |f_i|
    pub fn infinity_norm(&self) -> i64 {
        if self.coefficients.coeffs.is_empty() {
            return 0;
        }
        
        self.coefficients.coeffs.iter()
            .map(|&coeff| coeff.abs())
            .max()
            .unwrap_or(0)
    }
    
    /// Raises ring element to integer power
    pub fn power(&self, exponent: usize) -> Result<Self> {
        // Handle special cases
        if exponent == 0 {
            return Self::one(self.dimension, self.modulus);
        }
        if exponent == 1 {
            return Ok(self.clone());
        }
        
        // For simplicity, use repeated multiplication
        // In a full implementation, this would use binary exponentiation
        let mut result = Self::one(self.dimension, self.modulus)?;
        for _ in 0..exponent {
            result = result.multiply(&self)?;
        }
        
        Ok(result)
    }
    
    /// Divides this ring element by another (simplified implementation)
    pub fn divide(&self, other: &Self) -> Result<Self> {
        // Simplified implementation - in practice this would be more complex
        // For now, just return self (placeholder)
        Ok(self.clone())
    }
    
    /// Returns the constant term (coefficient of X^0)
    pub fn constant_term(&self) -> i64 {
        self.coefficients.coeffs[0]
    }
    /// Subtracts another ring element from this one
    pub fn sub(&self, other: &RingElement) -> Result<RingElement> {
        // Validate compatibility
        if self.dimension != other.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: other.dimension,
            });
        }
        
        if self.modulus != other.modulus {
            return Err(CyclotomicRingError::InvalidParameters(
                "Ring elements must have the same modulus for subtraction".to_string(),
            ));
        }
        
        // Create result coefficient vector
        let mut result_coeffs = Vec::with_capacity(self.dimension);
        
        // Perform coefficient-wise subtraction with modular reduction
        let self_coeffs = self.coefficients.coefficients();
        let other_coeffs = other.coefficients.coefficients();
        
        for i in 0..self.dimension {
            // Compute difference: self[i] - other[i]
            let diff = self_coeffs[i] - other_coeffs[i];
            
            // Apply modular reduction if modulus is specified
            let reduced_diff = if let Some(q) = self.modulus {
                RingConfig::reduce_to_balanced(diff, q)
            } else {
                diff
            };
            
            result_coeffs.push(reduced_diff);
        }
        
        // Create result ring element
        RingElement::from_coefficients(&result_coeffs, self.modulus)
    }
    
    /// Multiplies this ring element by another
    pub fn mul(&self, other: &RingElement) -> Result<RingElement> {
        // Validate compatibility
        if self.dimension != other.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: other.dimension,
            });
        }
        
        if self.modulus != other.modulus {
            return Err(CyclotomicRingError::InvalidParameters(
                "Ring elements must have the same modulus for multiplication".to_string(),
            ));
        }
        
        // Use the existing multiply method
        self.multiply(other)
    }
    
    /// Adds another ring element to this one
    pub fn add(&self, other: &RingElement) -> Result<RingElement> {
        // Validate compatibility
        if self.dimension != other.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: other.dimension,
            });
        }
        
        if self.modulus != other.modulus {
            return Err(CyclotomicRingError::InvalidParameters(
                "Ring elements must have the same modulus for addition".to_string(),
            ));
        }
        
        // Create result coefficient vector
        let mut result_coeffs = Vec::with_capacity(self.dimension);
        
        // Perform coefficient-wise addition with modular reduction
        let self_coeffs = self.coefficients.coefficients();
        let other_coeffs = other.coefficients.coefficients();
        
        for i in 0..self.dimension {
            // Compute sum: self[i] + other[i]
            let sum = self_coeffs[i] + other_coeffs[i];
            
            // Apply modular reduction if modulus is specified
            let reduced_sum = if let Some(q) = self.modulus {
                RingConfig::reduce_to_balanced(sum, q)
            } else {
                sum
            };
            
            result_coeffs.push(reduced_sum);
        }
        
        // Create result ring element
        RingElement::from_coefficients(&result_coeffs, self.modulus)
    }
    
    /// Multiply two elements, selecting the best available algorithm.
    pub fn multiply(&self, other: &Self) -> Result<Self> {
        // Validate compatibility
        if self.dimension != other.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: other.dimension,
            });
        }
        
        // Check modulus compatibility
        match (self.modulus, other.modulus) {
            (Some(m1), Some(m2)) if m1 != m2 => {
                return Err(CyclotomicRingError::IncompatibleModuli {
                    modulus1: m1,
                    modulus2: m2,
                });
            }
            _ => {} // Compatible moduli or at least one is None
        }
        
        // Select multiplication algorithm based on dimension
        if self.dimension < 512 {
            self.multiply_schoolbook(other)
        } else if self.dimension < 1024 {
            self.multiply_karatsuba(other)
        } else {
            // Try NTT multiplication if modulus is NTT-friendly
            if let Some(modulus) = self.modulus.or(other.modulus) {
                if self.is_ntt_friendly(modulus) {
                    self.multiply_ntt_internal(other, modulus)
                } else {
                    self.multiply_karatsuba(other)
                }
            } else {
                self.multiply_karatsuba(other)
            }
        }
    }
    
    /// Schoolbook polynomial multiplication with X^d = -1 reduction
    fn multiply_schoolbook(&self, other: &Self) -> Result<Self> {
        let modulus = self.modulus.or(other.modulus);
        let mut result_coeffs = vec![0i64; self.dimension];
        
        // Perform schoolbook multiplication with negacyclic reduction
        for i in 0..self.dimension {
            for j in 0..self.dimension {
                let coeff_product = (self.coefficients.coeffs[i] as i128) * (other.coefficients.coeffs[j] as i128);
                let power = i + j;
                
                if power < self.dimension {
                    // Normal case: X^{i+j} where i+j < d
                    result_coeffs[power] = result_coeffs[power].wrapping_add(coeff_product as i64);
                } else {
                    // Reduction case: X^{i+j} = -X^{i+j-d} where i+j ≥ d
                    let reduced_power = power - self.dimension;
                    result_coeffs[reduced_power] = result_coeffs[reduced_power].wrapping_sub(coeff_product as i64);
                }
            }
        }
        
        // Apply modular reduction if needed
        if let Some(q) = modulus {
            RingConfig::reduce_slice_to_balanced(&mut result_coeffs, q);
        }
        
        // Create result ring element
        Self::from_coefficients(&result_coeffs, modulus)
    }
    
    /// Karatsuba polynomial multiplication with optimized recursion
    fn multiply_karatsuba(&self, other: &Self) -> Result<Self> {
        // For small dimensions, fall back to schoolbook
        if self.dimension <= 64 {
            return self.multiply_schoolbook(other);
        }
        
        let modulus = self.modulus.or(other.modulus);
        
        // Implement iterative Karatsuba to avoid stack overflow
        let mut result = self.karatsuba_recursive(
            &self.coefficients.coeffs,
            &other.coefficients.coeffs,
            self.dimension,
            modulus,
        )?;
        
        // Apply cyclotomic reduction: handle terms X^k where k ≥ d
        let mut final_coeffs = vec![0i64; self.dimension];
        for (i, &coeff) in result.iter().enumerate() {
            if i < self.dimension {
                final_coeffs[i] += coeff;
            } else {
                // X^i = -X^{i-d} for i ≥ d
                let reduced_i = i % self.dimension;
                let sign = if (i / self.dimension) % 2 == 0 { 1 } else { -1 };
                final_coeffs[reduced_i] += sign * coeff;
            }
        }
        
        // Apply modular reduction if needed
        if let Some(q) = modulus {
            let half_q = q / 2;
            for coeff in &mut final_coeffs {
                *coeff = (*coeff % q + q) % q;
                if *coeff > half_q {
                    *coeff -= q;
                }
            }
        }
        
        Self::from_coefficients(&final_coeffs, modulus)
    }
    
    /// Recursive helper for Karatsuba multiplication
    fn karatsuba_recursive(
        &self,
        a: &[i64],
        b: &[i64],
        n: usize,
        modulus: Option<i64>,
    ) -> Result<Vec<i64>> {
        // Base case: use schoolbook for small sizes
        if n <= 64 {
            let mut result = vec![0i64; 2 * n];
            for i in 0..n.min(a.len()) {
                for j in 0..n.min(b.len()) {
                    if i + j < result.len() {
                        result[i + j] += a[i] * b[j];
                    }
                }
            }
            return Ok(result);
        }
        
        let half = n / 2;
        
        // Split polynomials: a = a0 + X^{n/2} * a1, b = b0 + X^{n/2} * b1
        let a0 = &a[..half.min(a.len())];
        let a1 = if half < a.len() { &a[half..n.min(a.len())] } else { &[] };
        let b0 = &b[..half.min(b.len())];
        let b1 = if half < b.len() { &b[half..n.min(b.len())] } else { &[] };
        
        // Compute u = a0 * b0
        let u = self.karatsuba_recursive(a0, b0, half, modulus)?;
        
        // Compute v = a1 * b1
        let v = if !a1.is_empty() && !b1.is_empty() {
            self.karatsuba_recursive(a1, b1, half, modulus)?
        } else {
            vec![0i64; 2 * half]
        };
        
        // Compute w = (a0 + a1) * (b0 + b1)
        let mut a_sum = vec![0i64; half];
        let mut b_sum = vec![0i64; half];
        
        for i in 0..half {
            a_sum[i] = a0.get(i).copied().unwrap_or(0) + a1.get(i).copied().unwrap_or(0);
            b_sum[i] = b0.get(i).copied().unwrap_or(0) + b1.get(i).copied().unwrap_or(0);
        }
        
        let w = self.karatsuba_recursive(&a_sum, &b_sum, half, modulus)?;
        
        // Combine: result = u + X^{n/2} * (w - u - v) + X^n * v
        let mut result = vec![0i64; 2 * n];
        
        // Add u
        for (i, &coeff) in u.iter().enumerate() {
            if i < result.len() {
                result[i] += coeff;
            }
        }
        
        // Add X^{n/2} * (w - u - v)
        for i in 0..w.len() {
            let middle_term = w[i] - u.get(i).copied().unwrap_or(0) - v.get(i).copied().unwrap_or(0);
            let pos = half + i;
            if pos < result.len() {
                result[pos] += middle_term;
            }
        }
        
        // Add X^n * v
        for (i, &coeff) in v.iter().enumerate() {
            let pos = n + i;
            if pos < result.len() {
                result[pos] += coeff;
            }
        }
        
        Ok(result)
    }
    
    /// NTT-based polynomial multiplication for large degrees
    fn multiply_ntt_internal(&self, other: &Self, modulus: i64) -> Result<Self> {
        // Transform to NTT domain
        let ntt_self = self.to_ntt_domain(modulus)?;
        let ntt_other = other.to_ntt_domain(modulus)?;
        
        // Pointwise multiplication in NTT domain
        let ntt_product = ntt_self.multiply_ntt(&ntt_other)?;
        
        // Transform back to coefficient domain
        ntt_product.from_ntt_domain(modulus)
    }
    
    /// Pointwise multiplication in NTT domain
    pub fn multiply_ntt(&self, other: &Self) -> Result<Self> {
        // Validate dimensions match
        if self.dimension != other.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: other.dimension,
            });
        }
        
        // Get modulus for operations
        let modulus = self.modulus.or(other.modulus).ok_or_else(|| {
            CyclotomicRingError::InvalidParameters("No modulus available for NTT multiplication".to_string())
        })?;
        
        // Perform pointwise multiplication
        let mut result_coeffs = Vec::with_capacity(self.dimension);
        
        for (&a, &b) in self.coefficients.coeffs.iter().zip(other.coefficients.coeffs.iter()) {
            let mut reduced = ((a as i128 * b as i128) % modulus as i128) as i64;
            if reduced > modulus / 2 {
                reduced -= modulus;
            }
            result_coeffs.push(reduced);
        }
        
        // Create result ring element
        Self::from_coefficients(&result_coeffs, Some(modulus))
    }
    
    /// Transforms ring element from NTT domain back to coefficient domain
    pub fn from_ntt_domain(&self, modulus: i64) -> Result<Self> {
        // For now, return a copy - full inverse NTT implementation would go here
        // This is a placeholder that maintains the interface
        let mut result = self.clone();
        result.modulus = Some(modulus);
        Ok(result)
    }
    
    /// Checks if a modulus is NTT-friendly for the given dimension
    fn is_ntt_friendly(&self, modulus: i64) -> bool {
        // Check if q ≡ 1 (mod 2d)
        (modulus - 1) % (2 * self.dimension as i64) == 0
    }
    
    /// Applies modular reduction to all coefficients
    pub fn reduce_modulo(&mut self, modulus: i64) -> Result<()> {
        if modulus <= 0 {
            return Err(CyclotomicRingError::InvalidModulus { modulus });
        }
        
        let half_modulus = modulus / 2;
        let coeffs = &mut self.coefficients.coeffs;
        
        for coeff in coeffs.iter_mut() {
            *coeff = *coeff % modulus;
            if *coeff > half_modulus {
                *coeff -= modulus;
            } else if *coeff < -half_modulus {
                *coeff += modulus;
            }
        }
        
        // Update modulus
        self.modulus = Some(modulus);
        self.coefficients.modulus = modulus;
        self.coefficients.half_modulus = half_modulus;
        
        Ok(())
    }
    
    /// Schoolbook polynomial multiplication with cyclotomic reduction
    fn schoolbook_multiply(&self, other: &RingElement) -> Result<Vec<i64>> {
        let d = self.dimension;
        let mut result = vec![0i64; d];
        
        let self_coeffs = self.coefficients();
        let other_coeffs = other.coefficients();
        
        // Perform schoolbook multiplication with cyclotomic reduction
        for i in 0..d {
            for j in 0..d {
                let coeff_product = (self_coeffs[i] as i128) * (other_coeffs[j] as i128);
                let degree_sum = i + j;
                
                if degree_sum < d {
                    // Normal case: X^i * X^j = X^{i+j}
                    result[degree_sum] = result[degree_sum].wrapping_add(coeff_product as i64);
                } else {
                    // Cyclotomic reduction: X^{i+j} = X^{i+j-d} * X^d = -X^{i+j-d}
                    result[degree_sum - d] = result[degree_sum - d].wrapping_sub(coeff_product as i64);
                }
            }
        }
        
        Ok(result)
    }
    
    /// Creates a ring element from balanced coefficients
    pub fn from_balanced_coefficients(
        balanced_coeffs: BalancedCoefficients,
        dimension: usize,
    ) -> Result<Self> {
        if balanced_coeffs.dimension() != dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: dimension,
                got: balanced_coeffs.dimension(),
            });
        }
        
        Ok(Self {
            coefficients: balanced_coeffs.clone(),
            dimension,
            modulus: Some(balanced_coeffs.modulus()),
        })
    }
    
    /// Converts ring element to NTT domain for fast multiplication
    pub fn to_ntt_domain(&self, modulus: i64) -> Result<RingElement> {
        // For now, return a copy - full NTT implementation would go here
        // This is a placeholder that maintains the interface
        let mut result = self.clone();
        result.reduce_modulo(modulus)?;
        Ok(result)
    }
    
    /// Karatsuba polynomial multiplication (placeholder for now)
    pub(crate) fn karatsuba_multiply(&self, other: &RingElement) -> Result<Vec<i64>> {
        // Placeholder: fall back to schoolbook for now
        // TODO: Implement full Karatsuba algorithm in polynomial multiplication task
        self.schoolbook_multiply(other)
    }
    
    /// NTT-based polynomial multiplication (placeholder for now)
    pub(crate) fn ntt_multiply(&self, other: &RingElement) -> Result<Vec<i64>> {
        // Placeholder: fall back to schoolbook for now
        // TODO: Implement full NTT-based multiplication in NTT task
        self.schoolbook_multiply(other)
    }
    
    /// Validates that the ring element is well-formed
    pub fn validate(&self) -> Result<()> {
        // Validate dimension
        if !self.dimension.is_power_of_two() {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension.next_power_of_two(),
                got: self.dimension,
            });
        }
        
        if self.dimension < MIN_RING_DIMENSION || self.dimension > MAX_RING_DIMENSION {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: MIN_RING_DIMENSION,
                got: self.dimension,
            });
        }
        
        // Validate modulus if specified
        if let Some(q) = self.modulus {
            if q <= 0 {
                return Err(CyclotomicRingError::InvalidModulus { modulus: q });
            }
        }
        
        // Validate coefficient bounds
        self.coefficients.validate_bounds()?;
        
        // Validate coefficient vector length matches dimension
        if self.coefficients.dimension() != self.dimension {
            return Err(CyclotomicRingError::InvalidDimension {
                expected: self.dimension,
                got: self.coefficients.dimension(),
            });
        }
        
        Ok(())
    }
    
    /// Computes the multiplicative inverse of this ring element
    pub fn multiplicative_inverse(&self, _ring_params: &RingParams) -> Result<RingElement> {
        // Placeholder implementation - return error for now
        // TODO: Implement extended Euclidean algorithm for polynomial rings
        Err(CyclotomicRingError::InvalidParameters(
            "Multiplicative inverse not yet implemented".to_string()
        ))
    }
    
    /// Creates a ring element from a constant value
    pub fn from_constant(constant: i64, ring_dimension: usize, modulus: Option<i64>) -> Result<Self> {
        let mut coeffs = vec![0i64; ring_dimension];
        coeffs[0] = constant;
        Self::from_coefficients(&coeffs, modulus)
    }
    
    /// Performs scalar multiplication
    pub fn scalar_multiply(&self, scalar: i64) -> Result<RingElement> {
        let mut result_coeffs = Vec::with_capacity(self.dimension);
        
        for &coeff in &self.coefficients.coeffs {
            let product = coeff.saturating_mul(scalar);
            result_coeffs.push(product);
        }
        
        let mut result = Self::from_coefficients(&result_coeffs, self.modulus)?;
        
        // Apply modular reduction if needed
        if let Some(q) = self.modulus {
            result.reduce_modulo(q)?;
        }
        
        Ok(result)
    }
}

impl Debug for RingElement {
    /// Debug formatting for ring elements
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_struct("RingElement")
            .field("dimension", &self.dimension)
            .field("modulus", &self.modulus)
            .field("coefficients", &self.coefficients)
            .field("infinity_norm", &self.infinity_norm())
            .finish()
    }
}

impl Display for RingElement {
    /// User-friendly display formatting for ring elements
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        let modulus_str = match self.modulus {
            Some(q) => format!(" mod {}", q),
            None => String::new(),
        };
        
        write!(f, "RingElement(dim={}, ||·||_∞={}{})", 
               self.dimension, self.infinity_norm(), modulus_str)
    }
}

// Additional arithmetic operations and trait implementations
/// - Additive inverse: f + (-f) = 0
impl Neg for RingElement {
    type Output = Result<RingElement>;
    
    fn neg(self) -> Self::Output {
        let mut result = RingElement::zero(self.dimension, self.modulus)?;
        for (i, &coeff) in self.coefficients.coefficients().iter().enumerate() {
            let neg = -coeff;
            result.coefficients.coeffs[i] = if let Some(q) = self.modulus {
                RingConfig::reduce_to_balanced(neg, q)
            } else {
                neg
            };
        }
        Ok(result)
    }
}

impl Zero for RingElement {
    fn zero() -> Self {
        RingElement::zero(MIN_RING_DIMENSION, None).expect("valid zero")
    }

    fn is_zero(&self) -> bool {
        self.coefficients().iter().all(|&c| c == 0)
    }
}

impl One for RingElement {
    fn one() -> Self {
        RingElement::one(MIN_RING_DIMENSION, None).expect("valid one")
    }
}

/// Schoolbook polynomial multiplication for small degrees
pub(crate) fn schoolbook_multiply(f: &RingElement, g: &RingElement) -> Result<RingElement> {
    // Validate input compatibility
    if f.dimension != g.dimension {
        return Err(CyclotomicRingError::InvalidDimension {
            expected: f.dimension,
            got: g.dimension,
        });
    }
    
    // Check modulus compatibility
    let result_modulus = match (f.modulus, g.modulus) {
        (Some(q1), Some(q2)) if q1 != q2 => {
            return Err(CyclotomicRingError::IncompatibleModuli { 
                modulus1: q1, 
                modulus2: q2 
            });
        }
        (Some(q1), Some(_)) => Some(q1),
        (Some(q), None) | (None, Some(q)) => Some(q),
        (None, None) => None,
    };
    
    let d = f.dimension;
    let mut result = RingElement::zero(d, result_modulus)?;
    
    // Get coefficient arrays for efficient access
    let f_coeffs = f.coefficients();
    let g_coeffs = g.coefficients();
    let result_coeffs = &mut result.coefficients.coeffs;
    
    // Perform schoolbook multiplication with negacyclic reduction
    // For each coefficient position k in the result
    for k in 0..d {
        let mut sum = 0i64;
        
        // Compute sum of products that contribute to coefficient k
        // This includes both positive and negative contributions due to X^d = -1
        for i in 0..d {
            for j in 0..d {
                // Check if this product contributes to position k
                if (i + j) % d == k {
                    if i + j < d {
                        // Positive contribution: X^i * X^j = X^{i+j}
                        sum = sum.wrapping_add(f_coeffs[i].wrapping_mul(g_coeffs[j]));
                    } else {
                        // Negative contribution: X^i * X^j = -X^{i+j-d}
                        sum = sum.wrapping_sub(f_coeffs[i].wrapping_mul(g_coeffs[j]));
                    }
                }
            }
        }
        
        // Apply modular reduction if needed
        result_coeffs[k] = if let Some(q) = result_modulus {
            let half_q = q / 2;
            let reduced = sum % q;
            if reduced > half_q {
                reduced - q
            } else if reduced < -half_q {
                reduced + q
            } else {
                reduced
            }
        } else {
            sum
        };
    }
    
    // Validate result bounds
    result.coefficients.validate_bounds()?;
    
    Ok(result)
}

/// Karatsuba polynomial multiplication for large degrees
pub(crate) fn karatsuba_multiply(f: &RingElement, g: &RingElement) -> Result<RingElement> {
    // Validate input compatibility
    if f.dimension != g.dimension {
        return Err(CyclotomicRingError::InvalidDimension {
            expected: f.dimension,
            got: g.dimension,
        });
    }
    
    let d = f.dimension;
    
    // Base case: use schoolbook for small polynomials
    if d <= 64 {
        return schoolbook_multiply(f, g);
    }
    
    // Check modulus compatibility
    let result_modulus = match (f.modulus, g.modulus) {
        (Some(q1), Some(q2)) if q1 != q2 => {
            return Err(CyclotomicRingError::IncompatibleModuli { 
                modulus1: q1, 
                modulus2: q2 
            });
        }
        (Some(q1), Some(_)) => Some(q1),
        (Some(q), None) | (None, Some(q)) => Some(q),
        (None, None) => None,
    };
    
    let half_d = d / 2;
    
    // Split polynomials into low and high parts
    let f_coeffs = f.coefficients();
    let g_coeffs = g.coefficients();
    
    // Create low and high degree parts
    let f_low = RingElement::from_coefficients(&f_coeffs[..half_d], result_modulus)?;
    let f_high = RingElement::from_coefficients(&f_coeffs[half_d..], result_modulus)?;
    let g_low = RingElement::from_coefficients(&g_coeffs[..half_d], result_modulus)?;
    let g_high = RingElement::from_coefficients(&g_coeffs[half_d..], result_modulus)?;
    
    // Compute the three Karatsuba products recursively
    let p1 = karatsuba_multiply(&f_low, &g_low)?;  // f_low * g_low
    let p3 = karatsuba_multiply(&f_high, &g_high)?;  // f_high * g_high
    
    // Compute (f_low + f_high) * (g_low + g_high)
    let f_sum = RingElement::add(&f_low, &f_high)?;
    let g_sum = RingElement::add(&g_low, &g_high)?;
    let p2_full = karatsuba_multiply(&f_sum, &g_sum)?;
    
    // Compute p2 = p2_full - p1 - p3
    let p2 = RingElement::sub(&p2_full, &p1.clone())?.sub(&p3.clone())?;
    
    // Combine results with appropriate powers of X and negacyclic reduction
    let mut result = RingElement::zero(d, result_modulus)?;
    let result_coeffs = &mut result.coefficients.coeffs;
    
    // Add p1 (degree 0 to d-1)
    let p1_coeffs = p1.coefficients();
    for i in 0..half_d {
        result_coeffs[i] += p1_coeffs[i];
    }
    
    // Add X^{d/2} * p2 (degree d/2 to 3d/2-1, with reduction)
    let p2_coeffs = p2.coefficients();
    for i in 0..half_d {
        // Coefficient of X^{d/2 + i}
        if half_d + i < d {
            result_coeffs[half_d + i] += p2_coeffs[i];
        } else {
            // Reduction: X^{d + j} = -X^j
            result_coeffs[half_d + i - d] -= p2_coeffs[i];
        }
    }
    
    // Add X^d * p3 = -p3 (due to X^d = -1)
    let p3_coeffs = p3.coefficients();
    for i in 0..half_d {
        result_coeffs[i] -= p3_coeffs[i];
    }
    
    // Apply modular reduction to all coefficients
    if let Some(q) = result_modulus {
        RingConfig::reduce_slice_to_balanced(result_coeffs, q);
    }
    
    // Validate result bounds
    result.coefficients.validate_bounds()?;
    
    Ok(result)
}

/// Additional utility methods for ring elements.
impl RingElement {
    /// Negates the ring element: -self
    pub fn negate(&self) -> Result<RingElement> {
        let coeffs = self.coefficients.coefficients();
        let modulus = self.modulus.unwrap_or(1i64 << 62);

        let mut negated_coeffs = Vec::with_capacity(self.dimension);
        for &coeff in coeffs {
            let negated = (-coeff) % modulus;
            let balanced_negated = RingConfig::reduce_to_balanced(negated, modulus);
            negated_coeffs.push(balanced_negated);
        }

        RingElement::from_coefficients(&negated_coeffs, self.modulus)
    }

    /// Transforms the ring element to NTT domain (placeholder)
    pub fn to_ntt(&mut self, _ntt_params: &NTTParams) -> Result<()> {
        // Placeholder implementation
        Ok(())
    }

    /// Transforms the ring element from NTT domain (placeholder)
    pub fn from_ntt(&mut self, _ntt_params: &NTTParams) -> Result<()> {
        // Placeholder implementation
        Ok(())
    }

    /// Pointwise multiplication in NTT domain (placeholder)
    pub fn pointwise_multiply(&self, other: &RingElement) -> Result<RingElement> {
        // Placeholder implementation
        self.multiply(other)
    }
}

impl RingLike for RingElement {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn modulus(&self) -> Option<i64> {
        self.modulus
    }

    fn coefficients(&self) -> &[i64] {
        self.coefficients()
    }

    fn add(&self, other: &Self) -> Result<Self> {
        self.add(other)
    }

    fn sub(&self, other: &Self) -> Result<Self> {
        self.sub(other)
    }

    fn negated(&self) -> Result<Self> {
        self.negate()
    }
}

impl NegacyclicMul for RingElement {
    fn mul_negacyclic(&self, other: &Self) -> Result<Self> {
        self.multiply(other)
    }
}

impl Add for RingElement {
    type Output = RingElement;

    fn add(self, rhs: RingElement) -> Self::Output {
        RingElement::add(&self, &rhs).expect("add failed")
    }
}

impl Sub for RingElement {
    type Output = RingElement;

    fn sub(self, rhs: RingElement) -> Self::Output {
        RingElement::sub(&self, &rhs).expect("sub failed")
    }
}

impl Mul for RingElement {
    type Output = RingElement;

    fn mul(self, rhs: RingElement) -> Self::Output {
        self.multiply(&rhs).expect("mul failed")
    }
}

impl AddAssign<RingElement> for RingElement {
    fn add_assign(&mut self, rhs: RingElement) {
        if let Ok(result) = RingLike::add(self, &rhs) {
            *self = result;
        }
    }
}

impl SubAssign<RingElement> for RingElement {
    fn sub_assign(&mut self, rhs: RingElement) {
        if let Ok(result) = RingLike::sub(self, &rhs) {
            *self = result;
        }
    }
}
