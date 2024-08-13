# Quantum Spin Chains
This crate is created by Emil Aagaard in 2024 as a part their Master's thesis "Quantum Spin Chains: Haldane's Conjecture and Symmetry-Protected Topological Phases".

The crate diagonalizes Hamiltonians on finite quantum spin chains with general S using power iteration. To decrease the dimensionality of the problem, it only considers the subspace of eigenvectors that
- have real coefficients in the natural basis.
- are eigenvectors of the total quantum spin operator in the z-direction.
- are time-reversal invariant.
- are reflection invariant.
- are translation invariant.
See the ``examples`` folder for calculation of the ground state energies and spectral gaps of the S=1 antiferromagnetic Heisenberg chain and the Affleck-Kennedy-Lieb-Tasaki chain.
