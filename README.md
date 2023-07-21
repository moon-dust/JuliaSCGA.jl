# JuliaSCGA
JuliaSCGA is a Julia package that implements the self-consistent Gaussian approximation (SCGA) method. Starting from a general spin Hamiltonian of isotropic or anisotropic exchange interactions and single ion anisotropy, the SCGA method calculates the spin-spin correlations in the designated reciprocal space. Compared to other methods like the classical Monte Carlo simulations, the SCGA method provides an efficient way to analyze the diffuse neutron scattering data, and can also be utilized to determine the spin Hamiltonian through fits to the experimental data.

In JuliaSCGA, the definition of the [UnitCell](src/UnitCell.jl), including the [InteractionMatrix](src/InteractionMatrix.jl) object, inherits from [SpinMC.jl](https://github.com/fbuessen/SpinMC.jl). As shown in the example, the J<sub>1</sub>-J<sub>2</sub> model on a diamond lattice can be defined as:

```julia
using JuliaSCGA

# define the basis vectors of the primitive cell
a1 = (7.5095, 0.0000, 0.0000) 
a2 = (3.7547, 6.5034, 0.0000) 
a3 = (3.7547, 2.1678, 6.1315) 
uc = UnitCell(a1,a2,a3) 

# add atom positions
b1 = addBasisSite!(uc, (1.8774, 1.0839, 0.7664)) 
b2 = addBasisSite!(uc, (13.1416, 7.5873, 5.3650)) 

# define coupling matrices and apply to the bonds
J1 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00]*-1.00 
J2 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00]* 0.25

addInteraction!(uc, b2, b1, J1, (1, 1, 0)) 
addInteraction!(uc, b2, b1, J1, (1, 0, 1)) 
addInteraction!(uc, b2, b1, J1, (0, 1, 1)) 
addInteraction!(uc, b2, b1, J1, (1, 1, 1)) 
addInteraction!(uc, b1, b1, J2, (1, -1, 0)) 
addInteraction!(uc, b1, b1, J2, (0, -1, 1)) 
addInteraction!(uc, b2, b2, J2, (1, -1, 0)) 
addInteraction!(uc, b2, b2, J2, (0, -1, 1)) 
addInteraction!(uc, b1, b1, J2, (1, 0, -1)) 
addInteraction!(uc, b2, b2, J2, (1, 0, -1)) 
addInteraction!(uc, b1, b1, J2, (1, 0, 0)) 
addInteraction!(uc, b1, b1, J2, (0, 1, 0)) 
addInteraction!(uc, b1, b1, J2, (0, 0, 1)) 
addInteraction!(uc, b2, b2, J2, (1, 0, 0)) 
addInteraction!(uc, b2, b2, J2, (0, 1, 0)) 
addInteraction!(uc, b2, b2, J2, (0, 0, 1)) 
```



