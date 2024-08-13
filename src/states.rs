//! This module contains the [`State`] struct,
//! which is a symmetrized version of a vector.
use crate::basis::BasisStates;
use rand::random;
use std::ops::{Mul, MulAssign, SubAssign};

/// Symmetrized version of a vector.
pub struct State {
    pub coefficients: Vec<f32>,
}

impl State {
    /// Creates the zero [`State`].
    pub fn from_zeros(basis_states_length: usize) -> Self {
        let coefficients = vec![0.0; basis_states_length];

        Self {coefficients}
    }

    /// Creates a random [`State`].
    pub fn from_random(basis_states_length: usize) -> Self {
        let coefficients = (0..basis_states_length)
            .map(|_| random())
            .collect();

        Self {coefficients}
    }

    /// Inner product between this [`State`] and another.
    pub fn dot(
        &self,
        state: &State,
        length: u8,
        basis_states: &BasisStates,
        symmetry_factors: &Vec<f32>
    ) -> f32 {
        let mut dot_product = 0.0;

        for (index, (coefficient, representer)) in self.coefficients.iter().zip(basis_states.representers.iter()).enumerate() {
            let symmetry_factor = symmetry_factors[representer.value];

            if symmetry_factor != 0.0 {
                dot_product += coefficient * state.coefficients[index] / symmetry_factor
            }
        }

        dot_product * 4.0 * length as f32
    }

    /// Gets the norm of the [`State`].
    fn get_norm(
        &self,
        length: u8,
        basis_states: &BasisStates,
        symmetry_factors: &Vec<f32>,
    ) -> f32 {
        self.dot(self, length, basis_states, symmetry_factors).sqrt()
    }

    /// Gets the infinity norm of the [`State`].
    pub fn get_infinity_norm(&self) -> f32 {
        self.coefficients.iter().max_by(|a, b| a.abs().total_cmp(&b.abs())).unwrap().abs()
    }

    /// Normalizes the [`State`].
    pub fn normalize(
        &mut self,
        length: u8,
        basis_states: &BasisStates,
        symmetry_factors: &Vec<f32>,
    ) {
        let norm = self.get_norm(length, basis_states, symmetry_factors);
        *self *= 1.0 / norm;
    }

    /// Normalizes the [`State`] using the infinity norm.
    pub fn infinity_normalize(&mut self) {
        let infinity_norm = self.get_infinity_norm();
        *self *= 1.0 / infinity_norm
    }

    /// Sets the [`State`] to zero.
    pub fn clear(&mut self) {
        for coefficient in self.coefficients.iter_mut() {
            *coefficient = 0.0;
        }
    }

    /// Converts the [`State`] into a regular vector.
    pub fn get_full_state(
        &self,
        basis_states: &BasisStates,
        symmetry_factors: &Vec<f32>
    ) -> Vec<f32> {
        let mut full_state = Vec::with_capacity(symmetry_factors.len());

        for (index, basis_state_index) in basis_states.representer_map.iter().enumerate() {
            let symmetry_factor = symmetry_factors[index];

            if symmetry_factor != 0.0 {
                full_state.push(self.coefficients[*basis_state_index])
            } else {
                full_state.push(0.0)
            }
        }
        
        full_state
    }
}

impl Mul<f32> for &State {
    type Output = State;
    /// Multiplies a [`State`] with a number.
    fn mul(self, scaler: f32) -> Self::Output {
        let mut coefficients = self.coefficients.clone();

        for coefficient in coefficients.iter_mut() {
            *coefficient *= scaler;
        }

        State {coefficients}
    }
}

impl MulAssign<f32> for State {
    /// Multiplies a [`State`] with a number.
    fn mul_assign(&mut self, scaler: f32) {
        for coefficient in self.coefficients.iter_mut() {
            *coefficient *= scaler;
        }
    }
}

impl SubAssign for State {
    /// Subtracts two [`State`]s.
    fn sub_assign(&mut self, state: Self) {
        for (index, coefficient) in self.coefficients.iter_mut().enumerate() {
            *coefficient -= state.coefficients[index];
        }
    }
}