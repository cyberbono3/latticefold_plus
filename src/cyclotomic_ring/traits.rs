use super::error::Result;

/// Minimal ring-like interface to decouple algorithms from storage.
pub trait RingLike: Sized {
    /// Returns the ring dimension.
    fn dimension(&self) -> usize;

    /// Returns the modulus if operating in Rq, otherwise None for integer ring.
    fn modulus(&self) -> Option<i64>;

    /// Borrow coefficients in balanced representation.
    fn coefficients(&self) -> &[i64];

    /// Add another element.
    fn add(&self, other: &Self) -> Result<Self>;

    /// Subtract another element.
    fn sub(&self, other: &Self) -> Result<Self>;

    /// Negate the element.
    fn negated(&self) -> Result<Self>;
}

/// Negacyclic multiplication interface for cyclotomic rings.
pub trait NegacyclicMul: Sized {
    /// Multiply two elements with X^d = -1 reduction.
    fn mul_negacyclic(&self, other: &Self) -> Result<Self>;
}
