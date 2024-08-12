use crate::basis::BasisStates;
use crate::states::State;
use crate::hamiltonians::{Hamiltonian, Implemented};

pub struct Model {
    pub base: u8,
    pub spin: f32,
    pub length: u8,
    pub total_s_z: u8,
    pub base_powers: Vec<usize>,
    pub flippers: Vec<isize>,
    pub basis_states: BasisStates,
    pub m_coefficients: Vec<f32>,
    pub p_coefficients: Vec<f32>,
}

impl Model {
    pub fn new(
        base: u8,
        length: u8,
        total_s_z: u8,
    ) -> Self {
        let spin = (base - 1) as f32 / 2.0;
        let base_powers = get_base_powers(base, length);
        let flippers = get_flippers(&base_powers);
        let basis_states = BasisStates::new(
            base,
            spin,
            length,
            total_s_z,
            &base_powers,
        );
        let m_coefficients = get_m_coefficients(base, spin);
        let p_coefficients = get_p_coefficients(base, spin);

        Self {
            base,
            spin,
            length,
            total_s_z,
            base_powers,
            flippers,
            basis_states,
            m_coefficients,
            p_coefficients,
        }
    }

    fn power_iterate<T>(
        &self,
        hamiltonian: &Hamiltonian<T>,
        iterations: u32,
        symmetry_factors: &Vec<f32>,
        lower_eigenpairs: &Vec<(State, f32)>,
    ) -> State where Hamiltonian<T>: Implemented {
        let mut state_0 = State::from_random(self.basis_states.length);
        let mut state_1 = State::from_zeros(self.basis_states.length);
        let max_eigenenergy = hamiltonian.get_max_eigenenergy(&self);

        for _ in 0..iterations {
            hamiltonian.apply(&state_0, &mut state_1, &self, symmetry_factors);
            state_1 -= &state_0 * max_eigenenergy;

            for (lower_eigenstate, lower_eigenenergy) in lower_eigenpairs.iter() {
                state_1 -= lower_eigenstate * (lower_eigenenergy * lower_eigenstate.dot(&state_0, self.length, &self.basis_states, &symmetry_factors));
            }

            state_0.clear();

            hamiltonian.apply(&state_1, &mut state_0, &self, symmetry_factors);
            state_0 -= &state_1 * max_eigenenergy;

            for (lower_eigenstate, lower_eigenenergy) in lower_eigenpairs.iter() {
                state_0 -= lower_eigenstate * (lower_eigenenergy * lower_eigenstate.dot(&state_1, self.length, &self.basis_states, &symmetry_factors));
            }

            state_1.clear();
            state_0.normalize(self.length, &self.basis_states, symmetry_factors);
        }
        
        state_0
    }

    pub fn find_eigenstate<T>(
        &mut self,
        hamiltonian: &Hamiltonian<T>,
        iterations: u32,
        symmetry_factors: &Vec<f32>,
        lower_eigenpairs: &Vec<(State, f32)>,
    ) -> (State, f32) where Hamiltonian<T>: Implemented {
        let eigenstate = self.power_iterate(hamiltonian, iterations, symmetry_factors, lower_eigenpairs);
        let mut eigenstate_times_eigenenergy = State::from_zeros(self.basis_states.length);
        hamiltonian.apply(&eigenstate, &mut eigenstate_times_eigenenergy, &self, symmetry_factors);
        let eigenenergy = eigenstate.dot(&eigenstate_times_eigenenergy, self.length, &self.basis_states, symmetry_factors);

        (eigenstate, eigenenergy)
    }
}

pub fn get_base_powers(base: u8, length: u8) -> Vec<usize> {
    (0..=length)
        .map(
            |power| (base as usize).pow((power % length) as u32)
        )
        .collect()
}

pub fn get_flippers(base_powers: &Vec<usize>) -> Vec<isize> {
    (0..base_powers.len() - 1)
        .map(
            |index| base_powers[index + 1] as isize - base_powers[index] as isize
        )
        .collect()
}

fn get_m_coefficients(base: u8, spin: f32) -> Vec<f32> {
    (0..base)
        .map(|index|
            (0.5 * (spin * (spin + 1.0) - (index as f32 - spin) * (index  as f32 - spin - 1.0)))
                .sqrt()
        )
        .collect()
}

fn get_p_coefficients(base: u8, spin: f32) -> Vec<f32> {
    (0..base)
        .map(|index|
            (0.5 * (spin * (spin + 1.0) - (index as f32 - spin) * (index  as f32 - spin + 1.0)))
                .sqrt()
        )
        .collect()
}