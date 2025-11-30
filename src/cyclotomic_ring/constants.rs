/// Extended Euclidean algorithm for computing gcd and Bézout coefficients.
#[allow(dead_code)]
pub(crate) fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (gcd, x1, y1) = extended_gcd(b, a % b);
        (gcd, y1, x1 - (a / b) * y1)
    }
}

/// Maximum supported ring dimension (power of 2).
pub const MAX_RING_DIMENSION: usize = 16384;

/// Minimum supported ring dimension (power of 2).
pub const MIN_RING_DIMENSION: usize = 1;

/// SIMD vector width for i64 operations (AVX-512 supports 8 x i64).
pub const SIMD_WIDTH: usize = 8;

/// Memory alignment for SIMD operations (64 bytes for AVX-512).
pub const MEMORY_ALIGNMENT: usize = 64;
