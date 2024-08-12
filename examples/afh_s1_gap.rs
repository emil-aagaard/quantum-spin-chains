use quantum_spin_chains::model::Model;
use quantum_spin_chains::hamiltonians::{Hamiltonian, Heisenberg};
use quantum_spin_chains::symmetries::get_symmetry_factors;

fn main() {
    // Defining parameters of the chain and the calculation
    let base = 3; // corresponds to S=1
    let length = 10; // chain length
    let total_s_z = 0; // searches only in the subspace with S^z_total=0
    let s = 1.0; // refers to s*H_AFH + (1-s)*H_trivial=H_AFH
    let iterations = 500; // iterations of power iteration
    
    // Defining the Hamiltonian and the model
    let hamiltonian = Hamiltonian::<Heisenberg>::new(s);
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

    // Defining storage for lower eigenvalues and corresponding 
    // eigenvectors (not relevant here)
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