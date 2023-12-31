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

# spin correlations in the k-plane defined in the cubic cell
k_grid = [0:0.005:3;]
qx = [x_grid for x_grid in k_grid, y_grid in k_grid]
qy = [y_grid for x_grid in k_grid, y_grid in k_grid]
q_calc = [vec(qx) vec(qy) zeros(length(qx))]

# transfer to the rlu of the primitive cell
qmat = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2]
q_crys = qmat*q_calc'
q_crys = q_crys'*1.0

# generate the bonding distance vector as a N*3 matrix
dist = getDist(uc)

# calculate the interaction matrix
Jq_calc = getFourier_iso(uc, dist, q_crys)

# solve lambda
β = 2.0
λ = solveLambda_iso(uc, β)

# finally, calculate the correlation function
correl = getCorr_iso(uc, Jq_calc, β, λ)

# plot diffuse pattern using GR for speed
using Plots
gr()
heatmap(qx[:,1], qy[1,:], reshape(correl, axes(qx)), xlabel="qx", ylabel="qy", aspect_ratio=1, size=(400,400))
savefig("diamond_0p25.png")


