use ark_ff::{AdditiveGroup, Fp, FpConfig};
use ark_std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_traits::Zero;

use super::traits::PolyRing;

/// Cyclotomic ring element backed by fixed-size coefficient array.
#[derive(Clone, PartialEq, Eq)]
pub struct CyclotomicRing<const N: usize, const D: usize, C: FpConfig<N>> {
    pub(crate) coefficients: [Fp<C, N>; D],
}

impl<const N: usize, const D: usize, C: FpConfig<N>> Add
    for CyclotomicRing<N, D, C>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut coeffs = self.coefficients;
        for (c, r) in coeffs.iter_mut().zip(rhs.coefficients.iter()) {
            *c += *r;
        }
        Self {
            coefficients: coeffs,
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> AddAssign
    for CyclotomicRing<N, D, C>
{
    fn add_assign(&mut self, rhs: Self) {
        for (c, r) in self.coefficients.iter_mut().zip(rhs.coefficients.iter())
        {
            *c += *r;
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> Sub
    for CyclotomicRing<N, D, C>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coeffs = self.coefficients;
        for (c, r) in coeffs.iter_mut().zip(rhs.coefficients.iter()) {
            *c -= *r;
        }
        Self {
            coefficients: coeffs,
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> SubAssign
    for CyclotomicRing<N, D, C>
{
    fn sub_assign(&mut self, rhs: Self) {
        for (c, r) in self.coefficients.iter_mut().zip(rhs.coefficients.iter())
        {
            *c -= *r;
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> Neg
    for CyclotomicRing<N, D, C>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut coeffs = self.coefficients;
        for c in coeffs.iter_mut() {
            *c = -*c;
        }
        Self {
            coefficients: coeffs,
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> Mul
    for CyclotomicRing<N, D, C>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // Negacyclic convolution: X^D = -1
        let mut result = [<Fp<C, N> as AdditiveGroup>::ZERO; D];
        for i in 0..D {
            for j in 0..D {
                let prod = self.coefficients[i] * rhs.coefficients[j];
                let idx = (i + j) % D;
                if i + j < D {
                    result[idx] += prod;
                } else {
                    result[idx] -= prod;
                }
            }
        }
        Self {
            coefficients: result,
        }
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> MulAssign
    for CyclotomicRing<N, D, C>
{
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = core::mem::take(self);
        let product = lhs * rhs;
        *self = product;
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> Zero
    for CyclotomicRing<N, D, C>
{
    fn zero() -> Self {
        Self {
            coefficients: [<Fp<C, N> as AdditiveGroup>::ZERO; D],
        }
    }

    fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }
}

// impl<const N: usize, const D: usize, C:FpConfig<N>> CyclotomicRing<N, D, C> {
//     pub(super) fn from_coeffs_vec(mut coeffs: Vec<Fp<C::BaseFieldConfig, N>>) -> Self {
//         // C::reduce_in_place(&mut coeffs);

//         // Self(coeffs.try_into().expect("Wrong length"))
//         todo!()
//     }

//     fn from_fn<F>(f: F) -> Self
//     // where
//     //     F: FnMut(usize) -> Fp<C::BaseFieldConfig, N>,
//     // {
//     //     let coeffs = core::array::from_fn::<_, D, _>(f);
//     //     Self::from_array(coeffs)
//     todo!()
//     }

//     pub(crate) fn from_array(coeffs: [Fp<C::BaseFieldConfig, N>; D]) -> Self {
//         Self(coeffs)
//     }

//     fn poly_mul(&self, rhs: &Self) -> Self {
//         let lhs_coeffs = self.coeffs();
//         let rhs_coeffs = rhs.coeffs();
//         let mut coeffs = vec![<Fp::<C::BaseFieldConfig, N> as Field>::ZERO; 2 * D - 1];
//         for (i, &lhs_coeff) in lhs_coeffs.iter().enumerate() {
//             if !lhs_coeff.is_zero() {
//                 for (j, &rhs_coeff) in rhs_coeffs.iter().enumerate() {
//                     coeffs[i + j] += lhs_coeff * rhs_coeff;
//                 }
//             }
//         }
//         C::reduce_in_place(&mut coeffs);
//         Self::from_coeffs_vec(coeffs)
//     }

//     fn poly_mul_unchecked(&self, rhs: &Self) -> Self {
//         let lhs_coeffs = self.coeffs();
//         let rhs_coeffs = rhs.coeffs();
//         let mut coeffs = vec![<Fp::<C::BaseFieldConfig, N> as Field>::ZERO; 2 * D - 1];
//         for (i, &lhs_coeff) in lhs_coeffs.iter().enumerate() {
//             for (j, &rhs_coeff) in rhs_coeffs.iter().enumerate() {
//                 coeffs[i + j] += lhs_coeff * rhs_coeff;
//             }
//         }
//         C::reduce_in_place(&mut coeffs);
//         Self::from_coeffs_vec(coeffs)
//     }

//     fn poly_mul_in_place(&mut self, rhs: &Self) {
//         // we need a reduce function for SVector to properly do a multiplication in place
//         let res = *self * rhs;
//         self.0 = res.0;
//     }

// }

impl<const N: usize, const D: usize, C: FpConfig<N>> Default
    for CyclotomicRing<N, D, C>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const N: usize, const D: usize, C: FpConfig<N>> PolyRing
    for CyclotomicRing<N, D, C>
{
    type BaseRing = Fp<C, N>;

    fn coeffs(&self) -> &[Self::BaseRing] {
        &self.coefficients
    }

    fn coeffs_mut(&mut self) -> &mut [Self::BaseRing] {
        &mut self.coefficients
    }

    fn into_coeffs(self) -> Vec<Self::BaseRing> {
        self.coefficients.to_vec()
    }

    fn dimension(&self) -> usize {
        D
    }
}
