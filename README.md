# Quantum Spin Chains
This crate is created by Emil Aagaard in 2024 as a part their Master's thesis "Quantum Spin Chains: Haldane's Conjecture and Symmetry-Protected Topological Phases".

The crate diagonalizes Hamiltonians on finite quantum spin chains with general <img src="https://latex.codecogs.com/svg.image?%5Cinline%20%5Csmall%20%7B%5Ccolor%7BWhite%7D%7DS%5Cin%5Cmathbb%20N%5Cslash2"/> using power iteration. To decrease the dimensionality of the problem, it only considers the subspace of eigenvectors that
- have real coefficients in the natural basis (i.e., the basis that such that <img src="https://latex.codecogs.com/svg.image?\bg{white}\hat&space;S_n^z\vert\dots\sigma_n\dots\rangle=\sigma_n\vert\dots\sigma_n\dots\rangle"/>).
- are eigenvectors of <img src="https://latex.codecogs.com/svg.image?\bg{white}\hat&space;S_{\text{total}}^z"/>.
- are time-reversal invariant.
- are reflection invariant.
- are translation invariant.
See the examples folder for calculation of the ground state energies and spectral gaps of the <img src="https://latex.codecogs.com/svg.image?\bg{white}S=1"/> antiferromagnetic Heisenberg chain and the Affleck-Kennedy-Lieb-Tasaki chain.