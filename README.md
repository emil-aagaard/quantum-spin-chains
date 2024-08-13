# Quantum Spin Chains
This crate is created by Emil Aagaard in 2024 as a part their Master's thesis "Quantum Spin Chains: Haldane's Conjecture and Symmetry-Protected Topological Phases".

The crate diagonalizes Hamiltonians on finite quantum spin chains with general spin using power iteration. It only considers the subspace of eigenvectors that
- have real coefficients in the natural basis (<img src="https://latex.codecogs.com/svg.image?\hat&space;S_n^z"/>).
- are 
- are time-reversal invariant.
- are reflection invariant.
- are translation invariant.