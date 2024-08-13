//! This Rust crate was created by Emil Aagaard in 2024, as part of the Master's thesis "Quantum Spin Chains: Haldane's Conjecture and Symmetry-Protected Topological Phases."
//! 
//! It diagonalizes Hamiltonians on finite quantum spin chains with general S using power iteration. To decrease the dimensionality of the problem, it only considers the subspace of eigenvectors that
//! - have real coefficients in the natural basis.
//! - are eigenvectors of the total quantum spin operator in the z-direction.
//! - are time-reversal invariant.
//! - are reflection invariant.
//! - are translation invariant.
//! 
//! See the [``examples``](https://github.com/emil-aagaard/quantum-spin-chains/tree/main/examples) folder for calculation of the ground state energies and spectral gaps of the S=1 antiferromagnetic Heisenberg (AFH) chain and the Affleck-Kennedy-Lieb-Tasaki (AKLT) chain.
pub mod symmetries;
pub mod basis;
pub mod states;
pub mod hamiltonians;
pub mod model;