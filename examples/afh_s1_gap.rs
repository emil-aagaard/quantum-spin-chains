//! In this example, the spectral gap of the S=1 AFH chain of length 10 
//! is calculated. It can be shown, that both the ground state and at 
//! least one first excited state is in the subspace where total 
//! quantum spin in the z-direction is 0. Also, it can be shown that 
//! the ground state is symmetric (has eigenvalue 1) under the three 
//! symmetries we consider, while the first excited state is 
//! antisymmetric (has eigenvalue -1).
use quantum_spin_chains::model::Model;
use quantum_spin_chains::hamiltonians::{Hamiltonian, AFH};
use quantum_spin_chains::symmetries::get_symmetry_factors;

fn main() {
    // Defining parameters of the chain and the calculation
    let base = 3; // corresponds to S=1
    let length = 10; // chain length
    let total_s_z = 0; // total quantum spin in the z-direction is 0
    let s = 1.0; // refers to s*H_AFH + (1-s)*H_triv=H_AFH
    let iterations = 500; // iterations of power iteration
    
    // Defining the Hamiltonian and the model
    let hamiltonian = Hamiltonian::<AFH>::new(s);
    let mut model = Model::new(base, length, total_s_z);

    // Symmetry eigenvalues for the ground state
    let symmetric_time_reversal_eigenvalue = false; // 1
    let symmetric_reflection_eigenvalue = false; // 1
    let symmetric_translation_eigenvalue = false; // 1
    let symmetric_symmetry_factors = get_symmetry_factors(
        &model.basis_states.symmetry_signs,
        symmetric_time_reversal_eigenvalue,
        symmetric_reflection_eigenvalue,
        symmetric_translation_eigenvalue
    );

    // Symmetry eigenvalues for the first excited state
    let antisymmetric_time_reversal_eigenvalue = true; // -1
    let antisymmetric_reflection_eigenvalue = true; // -1
    let antisymmetric_translation_eigenvalue = true; // -1
    let antisymmetric_symmetry_factors = get_symmetry_factors(
        &model.basis_states.symmetry_signs,
        antisymmetric_time_reversal_eigenvalue,
        antisymmetric_reflection_eigenvalue,
        antisymmetric_translation_eigenvalue
    );

    // Defining storage for lower eigenvalues and corresponding eigenvectors (not relevant here)
    let lower_eigenpairs = Vec::new();
    
    // Finding the ground state
    let (_gs, gs_energy) = model.find_eigenstate(
        &hamiltonian,
        iterations,
        &symmetric_symmetry_factors,
        &lower_eigenpairs
    );

    // Finding the first excited state
    let (_fes, fes_energy) = model.find_eigenstate(
        &hamiltonian,
        iterations,
        &antisymmetric_symmetry_factors,
        &lower_eigenpairs
    );

    println!("Ground state energy: {}", gs_energy);
    println!("First excited energy: {}", fes_energy);
    println!("Spectral gap: {}", fes_energy - gs_energy);
}