use std::fmt::{Display, Formatter, Result as FmtResult};

/// Result type for cyclotomic ring operations.
pub type Result<T> = std::result::Result<T, CyclotomicRingError>;

/// Errors that can occur in cyclotomic ring computations.
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum CyclotomicRingError {
    InvalidDimension { expected: usize, got: usize },
    InvalidModulus { modulus: i64 },
}

impl Display for CyclotomicRingError {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            CyclotomicRingError::InvalidDimension { expected, got } => {
                write!(
                    f,
                    "invalid dimension: expected {}, got {}",
                    expected, got
                )
            }
            CyclotomicRingError::InvalidModulus { modulus } => {
                write!(f, "invalid modulus: {}", modulus)
            }
        }
    }
}

impl std::error::Error for CyclotomicRingError {}
