using JuliaSCGA, LinearAlgebra
using Base.Threads, FFTW, Random

# define the square cell
a1 = (1.0000, 0.0000, 0.0000) 
a2 = (0.0000, 1.0000, 0.0000)
a3 = (0.0000, 0.0000, 1.0000)

uc = UnitCell(a1,a2,a3) 
b1 = addBasisSite!(uc, (0.0000, 0.0000, 0.000)) 

J1 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00]*-1.0
J2 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00]*0.5
J3 = [1.00 0.00 0.00; 0.00 1.00 0.00; 0.00 0.00 1.00]*0.25

addInteraction!(uc, b1, b1, J1, (1, 0, 0)) 
addInteraction!(uc, b1, b1, J1, (0, 1, 0)) 
addInteraction!(uc, b1, b1, J2, (1, -1, 0)) 
addInteraction!(uc, b1, b1, J2, (1, 1, 0)) 
addInteraction!(uc, b1, b1, J3, (2, 0, 0)) 
addInteraction!(uc, b1, b1, J3, (0, 2, 0)) 

# generate supercell and the k-grid
nH = Int(100)
nK = Int(100)
midH = Int(nH/2+1)
midK = Int(nK/2+1)

gridH = fftshift(fftfreq(nH))
gridK = fftshift(fftfreq(nK))

qmatH = [qh for qk in gridK, qh in gridH]
qmatK = [qk for qk in gridK, qh in gridH]
qlist = [qmatH[:] qmatK[:] zeros(length(qmatK[:]))]

dist = getDist(uc)
# calculate the interaction matrix
Jq = getFourier_iso(uc, dist, qlist)
Jq = real(dropdims(Jq, dims=(1,2)))
Jq = Jq .- minimum(Jq)
Jq_mat = reshape(Jq, nK, nH)

Ns = 3
V = nH*nK

# step1, define the initial values
Δ = 0.5
Σ_mat = -rand(Random.GLOBAL_RNG, nK, nH)*Jq_mat/10

# step2_1, calc Keff
Keff_mat = Jq_mat .+ Δ .- Σ_mat

T_new = 10
T_old = 1

while abs(T_new - T_old) > T_old*1e-6
# while T_new > -10

    global Keff_mat, Σ_mat, Δ, Keff_rM, Dq_rM, Dq_M_mP
    global T_new, T_old

    # step2_2, calc T
    T_old = T_new

    # calc Keff^-1 
    Keff_rM = Keff_mat.^-1

    # step3, calc Dq_rM
    Dq_rM = zeros(nK, nH)
    for indH = [1:nH;]
        for indK = [1:nK;]
            # order is [ -0.5, -0.4, -0.1, 0, 0.1, ..., 0.4,]
            Dq_rM[indK, indH] = Ns/2*sum(circshift(Keff_rM, (-indK+midK, -indH+midH)) .* Keff_rM)
        end
    end

    # step4, calc new Σ
    # flip the Dq_rM matrix for -p
    Dq_M_mP = reverse(Dq_rM.^-1, dims=(1,2))
    Dq_M_mP = circshift(Dq_M_mP, (1,1))

    # assing Dp at p = 0 to zero to exclude it from the summation
    Dq_M_mP[midK, midH] = 0

    Σ_mat = zeros(nK, nH)
    for indH = [1:nH;]
        for indK = [1:nK;]
            # order is [ -0.5, -0.4, -0.1, 0, 0.1, ..., 0.4,]
            Σ_mat[indK, indH] = -sum(circshift(Keff_rM, (-indK+midK, -indH+midH)) .* Dq_M_mP)
        end
    end

    Σmax = maximum(Σ_mat)
    Σ_mat = Σ_mat .- Σmax

    # update Keff and T
    Keff_mat = Jq_mat .+ Δ .- Σ_mat
    T_new = (Ns/(2*V)*sum(Keff_mat.^-1))^-1

    print(T_new,'\n')

end

using Plots
gr()
heatmap(gridH, gridK, Keff_rM, xlabel="qx", ylabel="qy", aspect_ratio=1, size=(400,400),xrange=(-0.52,0.5),yrange=(-0.52,0.5))
