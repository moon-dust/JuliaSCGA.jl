# JuliaSCGA
JuliaSCGA is a Julia package that implements the self-consistent Gaussian approximation (SCGA) method. Starting from a general spin Hamiltonian of isotropic or anisotropic exchange interactions and single ion anisotropy, the SCGA method calculates the spin-spin correlations in the designated reciprocal space. Compared to other methods like the classical Monte Carlo simulations, the SCGA method provides an efficient way to analyze the diffuse neutron scattering data, and can also be utilized to determine the spin Hamiltonian through fits to the experimental data.

The definition of the unit cell, including the [InteractionMatrix](src/InteractionMatrix.jl) object, inherits that of [SpinMC.jl](https://github.com/fbuessen/SpinMC.jl).



