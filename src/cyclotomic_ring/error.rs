use std::fmt::{Display, Formatter, Result as FmtResult};

/// Result type for cyclotomic ring operations.
pub type Result<T> = std::result::Result<T, CyclotomicRingError>;

/// Errors that can occur in cyclotomic ring computations.
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum CyclotomicRingError {
    InvalidDimension {
        expected: usize,
        got: usize,
    },
    InvalidModulus {
        modulus: i64,
    },
    CoefficientOutOfRange {
        coefficient: i64,
        min_bound: i64,
        max_bound: i64,
        position: usize,
    },
    InvalidParameters(String),
    IncompatibleModuli {
        modulus1: i64,
        modulus2: i64,
    },
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
            CyclotomicRingError::CoefficientOutOfRange {
                coefficient,
                min_bound,
                max_bound,
                position,
            } => write!(
                f,
                "coefficient {} out of range [{}, {}] at position {}",
                coefficient, min_bound, max_bound, position
            ),
            CyclotomicRingError::InvalidParameters(msg) => {
                write!(f, "invalid parameters: {}", msg)
            }
            CyclotomicRingError::IncompatibleModuli { modulus1, modulus2 } => {
                write!(f, "incompatible moduli: {} vs {}", modulus1, modulus2)
            }
        }
    }
}

impl std::error::Error for CyclotomicRingError {}
