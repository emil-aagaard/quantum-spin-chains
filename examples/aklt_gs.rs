//! In this example, the ground state energy of the AKLT 
//! chain of length 10 is calculated. It can be shown, 
//! that the ground state is in the subspace where total 
//! quantum spin in the z-direction is 0. Also, it can 
//! be shown that the is symmetric (has eigenvalue 1) 
//! under the three symmetries we consider.
use quantum_spin_chains::model::Model;
use quantum_spin_chains::hamiltonians::{Hamiltonian, AKLT};
use quantum_spin_chains::symmetries::get_symmetry_factors;

fn main() {
    // Defining parameters of the chain and the calculation
    let base = 3; // corresponds to S=1
    let length = 10; // chain length
    let total_s_z = 0; // total quantum spin in the z-direction is 0
    let time_reversal_eigenvalue = false; // 1
    let reflection_eigenvalue = false; // 1
    let translation_eigenvalue = false; // 1
    let s = 1.0; // refers to s*H_AKLT + (1-s)*H_triv=H_AKLT
    let iterations = 500; // iterations of power iteration
    
    // Defining the Hamiltonian and the model
    let hamiltonian = Hamiltonian::<AKLT>::new(s);
    let mut model = Model::new(base, length, total_s_z);
    let symmetry_factors = get_symmetry_factors(
        &model.basis_states.symmetry_signs,
        time_reversal_eigenvalue,
        reflection_eigenvalue,
        translation_eigenvalue
    );

    // Defining storage for lower eigenvalues and corresponding
    // eigenvectors (not relevant here)
    let lower_eigenpairs = Vec::new();
    
    // Finding the ground state
    let (_gs, gs_energy) = model.find_eigenstate(
        &hamiltonian,
        iterations,
        &symmetry_factors,
        &lower_eigenpairs
    );

    println!("Ground state energy: {}", gs_energy);
}