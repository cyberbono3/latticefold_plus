use num_bigint::BigUint;

pub trait MulUnchecked<Rhs = Self> {
    type Output;

    fn mul_unchecked(self, rhs: Rhs) -> Self::Output;
}

pub trait FromRandomBytes<T = Self>: Sized {
    fn byte_size() -> usize;

    fn try_from_random_bytes(bytes: &[u8]) -> Option<T>;
}

pub trait WithL2Norm {
    fn l2_norm_squared(&self) -> BigUint;
}

pub trait WithLinfNorm {
    fn linf_norm(&self) -> BigUint;
}
