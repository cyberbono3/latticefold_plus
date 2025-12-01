use ark_ff::{BitIteratorBE, Fp, FpConfig};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    One, UniformRand, Zero,
    fmt::{Debug, Display},
    hash::Hash,
    iter::{Product, Sum},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait Ring:
    'static
    + Copy
    + Clone
    + Debug
    + Display
    + Default
    + Send
    + Sync
    + Eq
    + Zero
    + One
    + Neg<Output = Self>
    + UniformRand
    + Sized
    + Hash
    + CanonicalSerialize
    + CanonicalDeserialize
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + Sum<Self>
    + Product<Self>
    + From<u128>
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
    + From<bool>
{
    /// The additive identity of the ring.
    const ZERO: Self;
    /// The multiplicative identity of the ring.
    const ONE: Self;

    /// Returns `sum([a_i * b_i])`.
    #[inline]
    fn sum_of_products<const T: usize>(a: &[Self; T], b: &[Self; T]) -> Self {
        let mut sum = Self::zero();
        for i in 0..a.len() {
            sum += a[i] * b[i];
        }
        sum
    }

    fn square_in_place(&mut self) -> &mut Self {
        *self *= *self;
        self
    }

    /// Returns `self^exp`, where `exp` is an integer represented with `u64` limbs,
    /// least significant limb first.
    #[must_use]
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();

        for i in BitIteratorBE::without_leading_zeros(exp) {
            res.square_in_place();

            if i {
                res *= *self;
            }
        }
        res
    }
}

impl<const N: usize, C: FpConfig<N>> Ring for Fp<C, N> {
    const ZERO: Self = <C as FpConfig<N>>::ZERO;
    const ONE: Self = <C as FpConfig<N>>::ONE;
}

pub trait PolyRing {
    type BaseRing: Ring;

    fn coeffs(&self) -> &[Self::BaseRing];
    fn coeffs_mut(&mut self) -> &mut [Self::BaseRing];
    fn into_coeffs(self) -> Vec<Self::BaseRing>;
    fn dimension(&self) -> usize;
}
