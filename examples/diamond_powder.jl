using JuliaSCGA, Random, Plots

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

# generate random wavevectors
function uniformOnSphere(rng = Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng) # rand(rng) generates a random number between 0 and 1
    z = 2.0 * rand(rng) - 1.0;
    r = sqrt(1.0 - z * z)
    return (r * cos(phi), r * sin(phi), z)
end

nRand::Integer = 5e3
Q_list = [0.1:0.1:5.5;]
Q_norm = reshape(transpose(repeat(Q_list, 1, nRand)), :,1)

q_lab = zeros(nRand*length(Q_list),3)
for ind_q = 1:size(q_lab,1)
    q_lab[ind_q,:] .= Q_norm[ind_q].*uniformOnSphere()
end

# transfer to the crystal coordinate system and calculate the interaction matrix
q_crys = q_lab*inv(uc.rl)
dist = getDist(uc)
Jq_calc = getFourier_iso(uc, dist, q_crys)

# solve lambda
temperature = 10 # unit in K
beta = 1/(temperature/11.6045)
lam_answer = solveLambda_iso(uc, beta)

# finally, calculate the correlation function and apply the magnetic form factor
correl = getCorr_iso(uc, Jq_calc, beta, lam_answer)

# magnetic form factor for Mn2+ spins
function FF_Mn2p(Q_vec::Matrix{Float64})
    # Q_vec is N-by-3 matrix of the q vector in the lab system in unit of \AA-1
    s_Fm = norm.(eachrow(Q_vec))/(4*pi)
    
    # <j0>
    A=0.4220; a=17.6840;
    B=0.5948; b=6.0050;
    C=0.0043; c=-0.6090;
    D=-0.0219;
    FF_J0 = A*exp.(-a*s_Fm.^2) .+ B*exp.(-b*s_Fm.^2) .+ C*exp.(-c*s_Fm.^2) .+ D;  
    FF2 = FF_J0.^2;

    return FF2
end

int_ff2 = correl .* FF_Mn2p(q_lab)
int_sum = dropdims(sum(reshape(int_ff2, nRand, :), dims=1), dims=1)./nRand

# plot diffuse pattern
gr()
Plots.plot(Q_list, int_sum, xlabel="Q (Å^-1)",ylabel="intensity (a.u.)",legend=false)
Plots.scatter!(Q_list, int_sum,legend=false)




