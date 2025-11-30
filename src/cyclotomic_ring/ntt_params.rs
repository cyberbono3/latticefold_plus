use super::RingConfig;
use super::error::{CyclotomicRingError, Result};

/// NTT parameters for transform-friendly moduli.
#[derive(Clone, Debug)]
pub struct NTTParams {
    pub dimension: usize,
    pub modulus: i64,
    pub root_of_unity: i64,
    pub twiddle_factors: Vec<i64>,
    pub bit_reversal_table: Vec<usize>,
}

impl NTTParams {
    /// Create parameters for power-of-two dimension with q ≡ 1 (mod 2d).
    pub fn new(dimension: usize, modulus: i64) -> Result<Self> {
        RingConfig::validate_dimension(dimension)?;
        if (modulus - 1) % (2 * dimension as i64) != 0 {
            return Err(CyclotomicRingError::InvalidModulus { modulus });
        }

        let root_of_unity = Self::find_primitive_root(modulus, 2 * dimension)?;

        let mut twiddle_factors = Vec::with_capacity(dimension);
        let mut power = 1i64;
        for _ in 0..dimension {
            twiddle_factors.push(power);
            power = (power * root_of_unity) % modulus;
        }

        let bit_reversal_table = Self::compute_bit_reversal_table(dimension);

        Ok(Self {
            dimension,
            modulus,
            root_of_unity,
            twiddle_factors,
            bit_reversal_table,
        })
    }

    fn find_primitive_root(modulus: i64, order: usize) -> Result<i64> {
        if (modulus - 1) % (order as i64) != 0 {
            return Err(CyclotomicRingError::InvalidParameters(format!(
                "Order {} does not divide p-1 = {}",
                order,
                modulus - 1
            )));
        }

        let generator = Self::find_generator(modulus)?;
        let exponent = (modulus - 1) / (order as i64);
        let root = Self::mod_pow(generator, exponent, modulus);

        if Self::mod_pow(root, order as i64, modulus) != 1 {
            return Err(CyclotomicRingError::InvalidParameters(
                "Failed to find primitive root of unity".to_string(),
            ));
        }

        Ok(root)
    }

    fn find_generator(p: i64) -> Result<i64> {
        for candidate in 2..p {
            if Self::is_generator(candidate, p) {
                return Ok(candidate);
            }
        }
        Err(CyclotomicRingError::InvalidParameters(format!(
            "No generator found for prime {}",
            p
        )))
    }

    fn is_generator(g: i64, p: i64) -> bool {
        let order = p - 1;
        if Self::mod_pow(g, order / 2, p) == 1 {
            return false;
        }
        true
    }

    fn mod_pow(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
        let mut result = 1;
        base %= modulus;
        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * base) % modulus;
            }
            exp >>= 1;
            base = (base * base) % modulus;
        }
        result
    }

    fn compute_bit_reversal_table(dimension: usize) -> Vec<usize> {
        let log_dim = dimension.trailing_zeros() as usize;
        let mut table = Vec::with_capacity(dimension);

        for i in 0..dimension {
            let mut reversed = 0;
            let mut temp = i;
            for _ in 0..log_dim {
                reversed = (reversed << 1) | (temp & 1);
                temp >>= 1;
            }
            table.push(reversed);
        }

        table
    }
}
