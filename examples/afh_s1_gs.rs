use quantum_spin_chains::model::Model;
use quantum_spin_chains::hamiltonians::{Hamiltonian, Heisenberg};
use quantum_spin_chains::symmetries::get_symmetry_factors;

fn main() {
    // Defining parameters of the chain and the calculation
    let base = 3; // corresponds to S=1
    let length = 10; // chain length
    let total_s_z = 0; // searches only in the subspace with S^z_total=0
    let time_reversal_eigenvalue = false; // 1
    let reflection_eigenvalue = false; // 1
    let translation_eigenvalue = false; // 1
    let s = 1.0; // refers to s*H_AFH + (1-s)*H_trivial=H_AFH
    let iterations = 500; // iterations of power iteration
    
    // Defining the Hamiltonian and the model
    let hamiltonian = Hamiltonian::<Heisenberg>::new(s);
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