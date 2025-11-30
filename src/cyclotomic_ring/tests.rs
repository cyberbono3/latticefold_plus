use super::*;
use crate::cyclotomic_ring::ring_element::{
    karatsuba_multiply, schoolbook_multiply,
};
use proptest::prelude::*;
use std::ops::Neg;

// Test data generators for property-based testing
prop_compose! {
    /// Generates valid ring dimensions (powers of 2 in supported range)
    fn valid_dimension()(exp in 5u32..15u32) -> usize {
        1 << exp  // 2^5 = 32 to 2^14 = 16384
    }
}

prop_compose! {
    /// Generates valid moduli for cryptographic applications
    fn valid_modulus()(q in 1000i64..1000000i64) -> i64 {
        // Ensure odd modulus for better mathematical properties
        if q % 2 == 0 { q + 1 } else { q }
    }
}

prop_compose! {
    /// Generates balanced coefficients within specified bounds
    fn balanced_coefficients(dimension: usize, modulus: i64)
                            (coeffs in prop::collection::vec(-modulus/2..=modulus/2, dimension)) -> Vec<i64> {
        coeffs
    }
}

#[test]
fn test_balanced_coefficients_creation() {
    // Test valid creation
    let coeffs = BalancedCoefficients::with_dimension(64, 1009).unwrap();
    assert_eq!(coeffs.dimension(), 64);
    assert_eq!(coeffs.modulus(), 1009);
    assert_eq!(coeffs.infinity_norm(), 0); // All coefficients are zero

    // Test invalid dimension (not power of 2)
    assert!(BalancedCoefficients::with_dimension(63, 1009).is_err());

    // Test invalid modulus (non-positive)
    assert!(BalancedCoefficients::with_dimension(64, 0).is_err());
    assert!(BalancedCoefficients::with_dimension(64, -1).is_err());
}

#[test]
fn test_standard_balanced_conversion() {
    let modulus = 1009i64;
    let standard_coeffs = vec![0, 1, 504, 505, 1008, 0, 0, 0]; // Power-of-two length

    // Convert to balanced representation
    let balanced =
        BalancedCoefficients::from_standard(&standard_coeffs, modulus).unwrap();

    // Expected balanced values: [0, 1, 504, -504, -1]
    let expected_balanced = vec![0, 1, 504, -504, -1, 0, 0, 0];
    assert_eq!(balanced.coefficients(), &expected_balanced);

    // Convert back to standard
    let recovered_standard = balanced.to_standard();
    assert_eq!(recovered_standard, standard_coeffs);
}

#[test]
fn test_ring_element_creation() {
    // Test zero element creation
    let zero = RingElement::zero(128, Some(1009)).unwrap();
    assert_eq!(zero.dimension(), 128);
    assert_eq!(zero.modulus(), Some(1009));
    assert_eq!(zero.constant_term(), 0);
    assert_eq!(zero.infinity_norm(), 0);

    // Test one element creation
    let one = RingElement::one(128, Some(1009)).unwrap();
    assert_eq!(one.constant_term(), 1);
    assert_eq!(one.infinity_norm(), 1);

    // Test creation from coefficients
    let coeffs = vec![1, 2, 3, 0, 0, 0, 0, 0]; // 8 coefficients
    let elem =
        RingElement::from_coefficients(&coeffs.clone(), Some(1009)).unwrap();
    assert_eq!(elem.dimension(), 8);
    assert_eq!(elem.constant_term(), 1);
}

proptest! {
    #[test]
    fn test_addition_properties(
        _dimension in valid_dimension(),
        modulus in valid_modulus(),
        coeffs1 in balanced_coefficients(32, 1009),  // Use fixed small dimension for speed
        coeffs2 in balanced_coefficients(32, 1009),
        coeffs3 in balanced_coefficients(32, 1009)
    ) {
        let f = RingElement::from_coefficients(&coeffs1, Some(modulus)).unwrap();
        let g = RingElement::from_coefficients(&coeffs2, Some(modulus)).unwrap();
        let h = RingElement::from_coefficients(&coeffs3, Some(modulus)).unwrap();
        let zero = RingElement::zero(32, Some(modulus)).unwrap();

        // Test commutativity: f + g = g + f
        let fg = f.clone().add(&g.clone()).unwrap();
        let gf = g.clone().add(&f.clone()).unwrap();
        prop_assert_eq!(fg.coefficients(), gf.coefficients());

        // Test associativity: (f + g) + h = f + (g + h)
        let fg_h = fg.add(&h.clone()).unwrap();
        let gh = g.add(&h).unwrap();
        let f_gh = f.clone().add(&gh).unwrap();
        prop_assert_eq!(fg_h.coefficients(), f_gh.coefficients());

        // Test identity: f + 0 = f
        let f_zero = f.clone().add(&zero).unwrap();
        prop_assert_eq!(f.coefficients(), f_zero.coefficients());

        // Test inverse: f + (-f) = 0
        let neg_f = f.clone().neg().unwrap();
        let f_neg_f = f.add(&neg_f).unwrap();
        prop_assert_eq!(f_neg_f.infinity_norm(), 0);
    }

    #[test]
    fn test_multiplication_properties(
        coeffs1 in balanced_coefficients(32, 1009),
        coeffs2 in balanced_coefficients(32, 1009),
        coeffs3 in balanced_coefficients(32, 1009)
    ) {
        let f = RingElement::from_coefficients(&coeffs1, Some(1009)).unwrap();
        let g = RingElement::from_coefficients(&coeffs2, Some(1009)).unwrap();
        let h = RingElement::from_coefficients(&coeffs3, Some(1009)).unwrap();
        let one = RingElement::one(32, Some(1009)).unwrap();

        // Test commutativity: f * g = g * f
        let fg = f.clone().mul(&g.clone()).unwrap();
        let gf = g.clone().mul(&f.clone()).unwrap();
        prop_assert_eq!(fg.coefficients(), gf.coefficients());

        // Test identity: f * 1 = f
        let f_one = f.clone().mul(&one).unwrap();
        prop_assert_eq!(f.coefficients(), f_one.coefficients());

        // Test distributivity: f * (g + h) = f * g + f * h
        let gh = g.clone().add(&h.clone()).unwrap();
        let f_gh = f.clone().mul(&gh).unwrap();
        let fg = f.clone().mul(&g).unwrap();
        let fh = f.mul(&h).unwrap();
        let fg_fh = fg.add(&fh).unwrap();
        prop_assert_eq!(f_gh.coefficients(), fg_fh.coefficients());
    }

    #[test]
    fn test_schoolbook_karatsuba_equivalence(
        coeffs1 in balanced_coefficients(64, 1009),
        coeffs2 in balanced_coefficients(64, 1009)
    ) {
        let f = RingElement::from_coefficients(&coeffs1, Some(1009)).unwrap();
        let g = RingElement::from_coefficients(&coeffs2, Some(1009)).unwrap();

        // Compute using both algorithms
        let schoolbook_result = schoolbook_multiply(&f, &g).unwrap();
        let karatsuba_result = karatsuba_multiply(&f, &g).unwrap();

        // Results should be identical
        prop_assert_eq!(schoolbook_result.coefficients(), karatsuba_result.coefficients());
    }
}

#[test]
fn test_negacyclic_property() {
    // Test that X^d = -1 in the ring
    let d = 8;
    let modulus = Some(1009);

    // Create X (polynomial with coefficient 1 for X^1)
    let mut x_coeffs = vec![0i64; d];
    x_coeffs[1] = 1; // X = 0 + 1*X + 0*X^2 + ...
    let x = RingElement::from_coefficients(&x_coeffs, modulus).unwrap();

    // Compute X^d by repeated multiplication
    let mut x_power = RingElement::one(d, modulus).unwrap();
    for _ in 0..d {
        x_power = x_power.mul(&x.clone()).unwrap();
    }

    // X^d should equal -1
    let neg_one = RingElement::one(d, modulus).unwrap().neg().unwrap();
    assert_eq!(x_power.coefficients(), neg_one.coefficients());
}

#[test]
fn test_simd_optimization() {
    // Test that SIMD and scalar implementations give same results
    let d = 64; // Multiple of SIMD_WIDTH
    let modulus = 1009i64;

    // Create test vectors with known values
    let coeffs1: Vec<i64> = (0..d as i64).map(|i| i % (modulus / 2)).collect();
    let coeffs2: Vec<i64> =
        (0..d as i64).map(|i| (i * 2) % (modulus / 2)).collect();

    let f = RingElement::from_coefficients(&coeffs1, Some(modulus)).unwrap();
    let g = RingElement::from_coefficients(&coeffs2, Some(modulus)).unwrap();

    // Test addition (uses SIMD internally)
    let sum = f.clone().add(&g.clone()).unwrap();

    // Verify result by manual computation
    for i in 0..d {
        let expected = (f.coefficients()[i] + g.coefficients()[i]) % modulus;
        let expected_balanced = if expected > modulus / 2 {
            expected - modulus
        } else {
            expected
        };
        assert_eq!(sum.coefficients()[i], expected_balanced);
    }
}

#[test]
fn test_memory_alignment() {
    // Test that coefficient vectors are properly aligned for SIMD
    let coeffs = BalancedCoefficients::with_dimension(128, 1009).unwrap();
    let ptr = coeffs.coefficients().as_ptr() as usize;

    // Check alignment (should be aligned to at least 8 bytes for i64)
    assert_eq!(ptr % 8, 0);

    // For optimal SIMD performance, should be aligned to larger boundaries
    // This is a soft requirement and may not always be achievable
}

#[test]
fn test_overflow_protection() {
    // Test behavior with large coefficients that might overflow
    let d = 32;
    let modulus = i64::MAX / 1000; // Large but not overflow-prone

    // Create coefficients near the modulus bound
    let large_coeffs: Vec<i64> = (0..d).map(|_| modulus / 2 - 1).collect();
    let f =
        RingElement::from_coefficients(&large_coeffs.clone(), Some(modulus))
            .unwrap();
    let g =
        RingElement::from_coefficients(&large_coeffs, Some(modulus)).unwrap();

    // Addition should not overflow and should maintain bounds
    let sum = f.add(&g).unwrap();
    assert!(sum.coefficients().iter().all(|&c| c.abs() <= modulus / 2));

    // Multiplication should also maintain bounds (though result may be reduced)
    let product = f.mul(&g).unwrap();
    assert!(
        product
            .coefficients()
            .iter()
            .all(|&c| c.abs() <= modulus / 2)
    );
}
