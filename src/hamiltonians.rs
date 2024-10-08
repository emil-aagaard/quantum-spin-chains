//! This module contains the objects needed to define 
//! a Hamiltonian. The Hamiltonians [`AFH`] and 
//! [`AKLT`] are already implemented.
use crate::states::State;
use crate::model::Model;

/// Defines the functionalities of a [`Hamiltonian`].
pub trait Implemented {
    /// Applying an operator.
    fn apply(
        &self,
        input_state: &State,
        output_state: &mut State,
        model: &Model,
        symmetry_factors: &Vec<f32>,
    );
    /// Gets the largest eigenvalue of an operator.
    fn get_max_eigenenergy(&self, model: &Model) -> f32;
}

/// Represents a Hamiltonian operator.
pub struct Hamiltonian<T> {
    /// The parameters of the Hamiltonian.
    pub parameters: T,
}

/// Holds the parameters a the AFH Hamiltonian.
pub struct AFH {
    /// s*H_AFH + (1-s)*H_triv.
    pub s: f32,
    /// 1-s.
    one_minus_s: f32,
}

/// Holds the parameters a the AKLT Hamiltonian.
pub struct AKLT {
    /// s*H_AKLT + (1-s)*H_triv.
    pub s: f32,
    /// 1-s.
    one_minus_s: f32,
}

impl Hamiltonian<AFH> {
    /// Creates the AFH Hamiltonian.
    pub fn new(s: f32) -> Self {
        let one_minus_s = 1.0 - s;

        Self {
            parameters: AFH {s, one_minus_s},
        }
    }
}

impl Hamiltonian<AKLT> {
    /// Creates the AKLT Hamiltonian.
    pub fn new(s: f32) -> Self {
        let one_minus_s = 1.0 - s;

        Self {
            parameters: AKLT {s, one_minus_s},
        }
    }
}

impl Implemented for Hamiltonian<AFH> {
    /// Applies the AFH Hamiltonian.
    fn apply(
        &self,
        input_state: &State,
        output_state: &mut State,
        model: &Model,
        symmetry_factors: &Vec<f32>
    ) {
        for basis_state_index in 0..model.basis_states.length {
            let representer = &model.basis_states.representers[basis_state_index];
            let symmetry_factor = symmetry_factors[representer.value];

            if symmetry_factor == 0.0 {
                continue;
            }

            let coefficient = input_state.coefficients[basis_state_index];
            let mut trivial_eigenvalue = 0.0;

            for chain_index in 0..model.length as usize {
                let digit = representer.digits[chain_index];
                let sigma = representer.sigmas[chain_index];
                let next_digit = representer.digits[chain_index + 1];
                let next_sigma = representer.sigmas[chain_index + 1];

                trivial_eigenvalue += sigma * (
                    next_sigma * self.parameters.s
                    + sigma * self.parameters.one_minus_s
                );

                if (digit != 0) && (next_digit != model.base - 1) {
                    let new_representer_value = (representer.value as isize + model.flippers[chain_index]) as usize;
                    let new_basis_state_index = model.basis_states.representer_map[new_representer_value];
                    let new_symmetry_factor = symmetry_factors[new_representer_value];
                    let symmetry_ratio = new_symmetry_factor / symmetry_factor;
                    let mp_coefficient = model.m_coefficients[digit as usize] * model.p_coefficients[next_digit as usize];

                    output_state.coefficients[new_basis_state_index] += coefficient 
                        * mp_coefficient
                        * symmetry_ratio
                        * self.parameters.s;
                }

                if (digit != model.base - 1) && (next_digit != 0) {
                    let new_representer_value = (representer.value as isize - model.flippers[chain_index]) as usize;
                    let new_basis_state_index = model.basis_states.representer_map[new_representer_value];
                    let new_symmetry_factor = symmetry_factors[new_representer_value];
                    let symmetry_ratio = new_symmetry_factor / symmetry_factor;
                    let pm_coefficient = model.p_coefficients[digit as usize] * model.m_coefficients[next_digit as usize];

                    output_state.coefficients[new_basis_state_index] += coefficient
                        * pm_coefficient
                        * symmetry_ratio
                        * self.parameters.s;
                }
            }

            output_state.coefficients[basis_state_index] += coefficient * trivial_eigenvalue
        }
    }

    fn get_max_eigenenergy(&self, model: &Model) -> f32 {
        model.spin.powi(2) * model.length as f32
    }
}

impl Implemented for Hamiltonian<AKLT> {
    /// Applies the AKLT Hamiltonian.
    fn apply(
        &self,
        input_state: &State,
        output_state: &mut State,
        model: &Model,
        symmetry_factors: &Vec<f32>
    ) {
        for basis_state_index in 0..model.basis_states.length {
            let representer = &model.basis_states.representers[basis_state_index];
            let symmetry_factor = symmetry_factors[representer.value];

            if symmetry_factor == 0.0 {
                continue;
            }

            let coefficient = input_state.coefficients[basis_state_index];
            let mut trivial_eigenvalue = 0.0;

            for chain_index in 0..model.length as usize {
                let flippers = project_2(
                    representer.digits[chain_index],
                    representer.digits[chain_index + 1],
                    model.flippers[chain_index],
                );

                trivial_eigenvalue += representer.sigmas[chain_index].powi(2);
                
                for (sign, flipper, cg_coefficient) in flippers {
                    if cg_coefficient != 0.0 {
                        let new_representer_value = if sign {
                            (representer.value as isize - flipper) as usize
                        } else {
                            (representer.value as isize + flipper) as usize
                        };

                        let new_basis_state_index = model.basis_states.representer_map[new_representer_value];
                        let new_symmetry_factor = symmetry_factors[new_representer_value];
                        let symmetry_ratio = new_symmetry_factor / symmetry_factor;

                        output_state.coefficients[new_basis_state_index] += coefficient
                            * cg_coefficient
                            * symmetry_ratio
                            * self.parameters.s;
                    } else {
                        break;
                    }
                }
            }

            output_state.coefficients[basis_state_index] += coefficient 
                * trivial_eigenvalue
                * self.parameters.one_minus_s;
        }
    }

    fn get_max_eigenenergy(&self, model: &Model) -> f32 {
        model.length as f32
    }
}

/// Projects onto the subspace with S=2.
fn project_2(
    digit: u8,
    next_digit: u8,
    flipper: isize,
) -> [(bool, isize, f32); 3] {
    match (digit, next_digit) {
        (0, 0) => [  // |--⟩
            (false, 0, 1.0), // 1*|--⟩
            (false, 0, 0.0),
            (false, 0, 0.0)
        ],
        (1, 0) => [ // |0-⟩
            (false, 0, 0.5), // 1/2*|0-⟩
            (false, flipper, 0.5), // 1/2*|-0⟩
            (false, 0, 0.0)
        ],
        (0, 1) => [ // |-0⟩
            (false, 0, 0.5), // 0.5*|-0⟩
            (true, flipper, 0.5), // 1/2*|0-⟩
            (false, 0, 0.0)
        ],
        (2, 0) => [ // |+-⟩
            (false, 0, 0.166666666), // 1/6*|+-⟩
            (false, flipper, 0.33333333), // 1/3*|00⟩
            (false, 2*flipper, 0.166666666) // 1/6*|-+⟩
        ],
        (1, 1) => [ // |00⟩
            (false, 0, 0.6666666), // 1/6*|00⟩
            (true, flipper, 0.33333333), // 1/3*|+-⟩
            (false, flipper, 0.33333333) // 1/3*|-+⟩
        ],
        (0, 2) => [ // |-+⟩
            (false, 0, 0.166666666), // 1/6*|-+⟩
            (true, flipper, 0.33333333), // 1/3*|00⟩
            (true, 2*flipper, 0.166666666) // 1/6*|+-⟩
        ],
        (2, 1) => [ // |+0⟩
            (false, 0, 0.5), // 1/2*|+0⟩
            (false, flipper, 0.5), // 1/2*|0+⟩
            (false, 0, 0.0)
        ],
        (1, 2) => [ // |0+⟩
            (false, 0, 0.5), // 1/2*|0+⟩
            (true, flipper, 0.5), // 1/2*|+0⟩
            (false, 0, 0.0)
        ],
        (2, 2) => [  // |++⟩
            (false, 0, 1.0), // 1*|++⟩
            (false, 0, 0.0),
            (false, 0, 0.0)
        ],
        _ => [
            (false, 0, 0.0),
            (false, 0, 0.0),
            (false, 0, 0.0)
        ]
    }
}


