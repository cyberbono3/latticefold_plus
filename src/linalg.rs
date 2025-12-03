use ark_std::vec::Vec;
use core::ops::{Add, Div, Index, Rem, Sub};
use num_integer::Integer;
use num_traits::Signed;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub type RowVector<T> = Vec<T>;
pub type Vector<T> = Vec<T>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<T>(pub Vec<RowVector<T>>);

impl<T> Matrix<T> {
    pub fn from_rows(rows: &[RowVector<T>]) -> Self
    where
        T: Clone,
    {
        Self(rows.to_vec())
    }

    pub fn row_iter(&self) -> core::slice::Iter<'_, RowVector<T>> {
        self.0.iter()
    }

    #[cfg(feature = "parallel")]
    pub fn par_row_iter(&self) -> rayon::slice::Iter<'_, RowVector<T>> {
        self.0.par_iter()
    }

    pub fn nrows(&self) -> usize {
        self.0.len()
    }

    pub fn ncols(&self) -> usize {
        self.0.first().map(|r| r.len()).unwrap_or(0)
    }
}

impl<T> From<Vec<RowVector<T>>> for Matrix<T> {
    fn from(rows: Vec<RowVector<T>>) -> Self {
        Self(rows)
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.0[index.0][index.1]
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymmetricMatrix<T>(pub Vec<Vec<T>>);

impl<T> SymmetricMatrix<T> {
    pub fn size(&self) -> usize {
        self.0.len()
    }
}

impl<T> From<Vec<Vec<T>>> for SymmetricMatrix<T> {
    fn from(rows: Vec<Vec<T>>) -> Self {
        Self(rows)
    }
}

impl<T> Index<(usize, usize)> for SymmetricMatrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (row, col) = if index.0 >= index.1 {
            (index.0, index.1)
        } else {
            (index.1, index.0)
        };

        &self.0[row][col]
    }
}

pub trait Transpose<T> {
    fn transpose(self) -> Vec<Vec<T>>;
}

impl<T: Clone> Transpose<T> for Vec<Vec<T>> {
    fn transpose(self) -> Vec<Vec<T>> {
        if self.is_empty() {
            return Vec::new();
        }

        let row_len = self[0].len();
        let mut transposed = vec![Vec::with_capacity(self.len()); row_len];

        for row in self {
            for (col_idx, value) in row.into_iter().enumerate() {
                if col_idx >= transposed.len() {
                    transposed.push(Vec::new());
                }
                transposed[col_idx].push(value);
            }
        }

        transposed
    }
}

pub fn rounded_div<T>(value: T, divisor: i128) -> T
where
    T: Integer
        + Signed
        + From<i128>
        + Add<i128, Output = T>
        + Sub<i128, Output = T>
        + Div<i128, Output = T>
        + Rem<i128, Output = T>
        + Clone,
{
    let half = T::from(divisor / 2);
    let quotient = value.clone() / divisor;
    let remainder = value % divisor;

    if remainder > half {
        quotient + 1
    } else if remainder < -half {
        quotient - 1
    } else {
        quotient
    }
}
