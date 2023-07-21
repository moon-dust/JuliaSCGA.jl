using JuliaSCGA, LinearAlgebra, Optim
using Base.Threads

# ShG230406, prototype SCGA code, J1-J2 model on the diamond lattice

# step1, define the model
a1 = (7.5095, 0.0000, 0.0000) 
a2 = (3.7547, 6.5034, 0.0000) 
a3 = (3.7547, 2.1678, 6.1315) 

uc = UnitCell(a1,a2,a3) 
b1 = addBasisSite!(uc, (1.8774, 1.0839, 0.7664)) 
b2 = addBasisSite!(uc, (13.1416, 7.5873, 5.3650)) 

J1 = [-1.00 -0.00 -0.00; -0.00 -1.00 -0.00; -0.00 -0.00 -1.00] 
J2 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00].*0.15

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


# define the reciprocal lattice, vectors are stored in rows
rl = 2*pi*inv([[a1...] [a2...] [a3...]])

# generate k_grid in the 1st BZ
n_grid = 8
k0 = [(1-n_grid)/(2*n_grid):1/n_grid:(n_grid-1)/(2*n_grid);]
qx = [x_grid for x_grid in k0, y_grid in k0, z_grid in k0]
qy = [y_grid for x_grid in k0, y_grid in k0, z_grid in k0]
qz = [z_grid for x_grid in k0, y_grid in k0, z_grid in k0]
q_BZ = [vec(qx) vec(qy) vec(qz)]

# interaction matrix in the 1st BZ
natom = length(uc.basis)
nbond = length(uc.interactions)
Jq = zeros(natom, natom, axes(q_BZ,1)).*im
for ind_q in axes(q_BZ,1)
    q_lab = q_BZ[ind_q,1].*rl[1,:] + q_BZ[ind_q,2].*rl[2,:] + q_BZ[ind_q,3].*rl[3,:]
    for ind_bond = 1:nbond
        local bond = uc.interactions[ind_bond]
        local pos1 = uc.basis[bond[1]]
        local pos2 = uc.basis[bond[2]] .+ bond[3][1].*uc.primitive[1] .+ bond[3][2].*uc.primitive[2] .+ bond[3][3].*uc.primitive[3]
        local dist = pos2 .- pos1
        Jq[bond[1], bond[2], ind_q] += bond[4][1]/2*exp(1im*(q_lab[1]*dist[1] + q_lab[2]*dist[2] + q_lab[3]*dist[3]))
        Jq[bond[2], bond[1], ind_q] += bond[4][1]/2*exp(-1im*(q_lab[1]*dist[1] + q_lab[2]*dist[2] + q_lab[3]*dist[3]))
    end
end

# eigensolution 
eig_val = zeros(natom, size(q_BZ,1))
for ind_q = 1:size(q_BZ,1)
    eig_val[:,ind_q] = eigvals(Jq[:,:,ind_q])
end

# solve lambda
beta = 8
lam_sum_chi2(lambda) = (sum(1 ./ (lambda[1] .+ beta.*eig_val))/size(eig_val,1)/size(q_BZ,1)-1/3)^2
lam_lb = maximum(-beta.*eig_val)

lower = lam_lb
upper = Inf
x0 = [lam_lb+10]
result = optimize(lam_sum_chi2, lower, upper, x0)
lam_answer = Optim.minimizer(result)[1]

# spin correlations in the k-plane defined in the cubic cell
k_grid = [0:0.02:3;]
qx = [x_grid for x_grid in k_grid, y_grid in k_grid]
qy = [y_grid for x_grid in k_grid, y_grid in k_grid]
q_calc = [vec(qx) vec(qy) zeros(length(qx))]

# transfer to the rlu of the primitive cell
qmat = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2]
q_crys = qmat*transpose(q_calc)
q_crys = transpose(q_crys)

# interaction matrix in the 1st BZ
Jq = zeros(natom, natom, axes(q_crys,1)).*im
Threads.@threads for ind_q in axes(q_crys,1)
    q_lab = q_crys[ind_q,1].*rl[1,:] + q_crys[ind_q,2].*rl[2,:] + q_crys[ind_q,3].*rl[3,:]
    for ind_bond = 1:nbond
        local bond = uc.interactions[ind_bond]
        local pos1 = uc.basis[bond[1]]
        local pos2 = uc.basis[bond[2]] .+ bond[3][1].*uc.primitive[1] .+ bond[3][2].*uc.primitive[2] .+ bond[3][3].*uc.primitive[3]
        local dist = pos2 .- pos1
        Jq[bond[1], bond[2], ind_q] += bond[4][1]/2*exp(1im*(q_lab[1]*dist[1] + q_lab[2]*dist[2] + q_lab[3]*dist[3]))
        Jq[bond[2], bond[1], ind_q] += bond[4][1]/2*exp(-1im*(q_lab[1]*dist[1] + q_lab[2]*dist[2] + q_lab[3]*dist[3]))
    end
end

# eigensolution 
# eig_val = zeros(natom, size(q_crys,1))
# for ind_q = 1:size(q_crys,1)
#     eig_val[:,ind_q] = eigvals(Jq[:,:,ind_q])
# end

correl = zeros(1, axes(q_crys,1)).*im
for ind_q in axes(q_crys,1)
    correl[ind_q] = sum(inv(lam_answer.*Diagonal(dropdims(ones(1,natom); dims=1)) .+ beta.*Jq[:,:,ind_q]));
end
correl = real(correl)

# plot diffuse pattern
using PyPlot
fig = PyPlot.figure(figsize=(6,6))
ax = fig.add_subplot()
pcontour = ax.pcolor(qx, qy, reshape(correl, axes(qx)))
xlabel("qx")
ylabel("qy")
# PyPlot.xlim(-1, 1)
# PyPlot.ylim(-1, 1)
title("diffuse pattern")
# ax.axes.grid(visible=false)
ax.set_aspect("equal")
colorbar(pcontour, fraction=0.046, pad=0.04)
display(gcf())



