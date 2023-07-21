using LinearAlgebra, Optim

struct UnitCell{D}
    primitive::NTuple{D,NTuple{D,Float64}}
    rl::Matrix{Float64}
    basis::Vector{NTuple{D,Float64}}
    interactions::Vector{Tuple{Int,Int,NTuple{D,Int},Matrix{Float64}}} #interactions specified as (basis1,basis2,offsetPrimitive,M)
    interactionsOnsite::Vector{Matrix{Float64}}
    interactionsField::Vector{Vector{Float64}}
    

    UnitCell(a1::NTuple{1,Float64}) = new{1}((a1,), Vector{NTuple{1,Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{1,Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
    UnitCell(a1::NTuple{2,Float64}, a2::NTuple{2,Float64}) = new{2}((a1,a2), Vector{NTuple{2,Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{2,Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
    UnitCell(a1::NTuple{3,Float64}, a2::NTuple{3,Float64}, a3::NTuple{3,Float64}) = new{3}((a1,a2,a3), 
                                                                                           2*pi*inv([[a1...] [a2...] [a3...]]),     
                                                                                           Vector{NTuple{3,Float64}}(undef,0), 
                                                                                           Vector{Tuple{Int,Int,NTuple{3,Int},Matrix{Float64}}}(undef,0), 
                                                                                           Vector{Matrix{Float64}}(undef,0), 
                                                                                           Vector{Vector{Float64}}(undef,0))
    UnitCell(primitives...) = new{length(primitives)}(primitives, Vector{NTuple{length(primitives),Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{length(primitives),Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
end

"""
Adds an interaction between spin1 located at basis site `b1` of the given `unitcell` and spin2 at basis site `b2` in a unit cell that is offset by `offset` lattice vectors. 
The exchange energy is calculated as spin1'.M.spin2. 
"""
function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::Matrix{Float64}, offset::NTuple{D,Int}=Tuple(zeros(Int,D))) where D
    size(M) == (3,3) || error("Interaction matrix must be of size 3x3.")
    b1 == b2 && offset == Tuple(zeros(Int,D)) && error("Interaction cannot be local. Use setInteractionOnsite!() instead.")

    push!(unitcell.interactions, (b1,b2,offset,M))
end

function setInteractionOnsite!(unitcell::UnitCell{D}, b::Int, M::Matrix{Float64}) where D
    size(M) == (3,3) || error("Interaction matrix must be of size 3x3.")
    unitcell.interactionsOnsite[b] = M
end

function setField!(unitcell::UnitCell{D}, b::Int, B::Vector{Float64}) where D
    size(B) == (3,) || error("Field must be a vector of length 3.")
    unitcell.interactionsField[b] = B
end

function addBasisSite!(unitcell::UnitCell{D}, position::NTuple{D,Float64}) where D
    push!(unitcell.basis, position)
    push!(unitcell.interactionsOnsite, zeros(3,3))
    push!(unitcell.interactionsField, zeros(3))
    return length(unitcell.basis)
end

function getDist(unitcell::UnitCell{D}) where D
    nbond = length(unitcell.interactions)
    dist = zeros(nbond,3)
    for ind_bond = 1:nbond
        bond = unitcell.interactions[ind_bond]
        pos1 = unitcell.basis[bond[1]]
        pos2 = unitcell.basis[bond[2]] .+ bond[3][1].*unitcell.primitive[1] .+ bond[3][2].*unitcell.primitive[2] .+ bond[3][3].*unitcell.primitive[3]
        dist[ind_bond,1] = pos2[1] .- pos1[1]
        dist[ind_bond,2] = pos2[2] .- pos1[2]
        dist[ind_bond,3] = pos2[3] .- pos1[3]
    end
    return dist
end

function getFourier_iso(unitcell::UnitCell{D}, dist::Matrix{Float64}, q_crys::Matrix{Float64}) where D
    natom = length(unitcell.basis)
    Jq = zeros(natom, natom, axes(q_crys,1)).*im
    Threads.@threads for ind_q in axes(q_crys,1)
        q_lab = zeros(3)
        q_lab[1] = q_crys[ind_q,1]*unitcell.rl[1,1] + q_crys[ind_q,2]*unitcell.rl[2,1] + q_crys[ind_q,3]*unitcell.rl[3,1]
        q_lab[2] = q_crys[ind_q,1]*unitcell.rl[1,2] + q_crys[ind_q,2]*unitcell.rl[2,2] + q_crys[ind_q,3]*unitcell.rl[3,2]
        q_lab[3] = q_crys[ind_q,1]*unitcell.rl[1,3] + q_crys[ind_q,2]*unitcell.rl[2,3] + q_crys[ind_q,3]*unitcell.rl[3,3]

        nbond = length(unitcell.interactions)
        for ind_bond = 1:nbond
            bond = unitcell.interactions[ind_bond]
            phase = q_lab[1]*dist[ind_bond, 1] + q_lab[2]*dist[ind_bond,2] + q_lab[3]*dist[ind_bond,3]
            Jterm =  bond[4][1]*exp(1im*phase) # no 1/2 factor, check logbook SCGA/A_further_check
            Jq[bond[1], bond[2], ind_q] += Jterm
            Jq[bond[2], bond[1], ind_q] += conj(Jterm)
        end
    end
    return Jq
end

function solveLambda_iso(unitcell::UnitCell{D}, beta::Float64) where D
    n_grid::Int64 = 8
    k0::Vector{Float64} = [(1-n_grid)/(2*n_grid):1/n_grid:(n_grid-1)/(2*n_grid);]
    qx::Array{Float64, 3} = [x_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    qy::Array{Float64, 3} = [y_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    qz::Array{Float64, 3} = [z_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    q_BZ = [vec(qx) vec(qy) vec(qz)]

    dist = getDist(unitcell)
    Jq_BZ = getFourier_iso(unitcell, dist, q_BZ)

    # eigensolution 
    natom = length(unitcell.basis)
    eig_val = zeros(natom, size(q_BZ,1))
    for ind_q in axes(q_BZ,1)
        eig_val[:,ind_q] = eigvals(Jq_BZ[:,:,ind_q])
    end

    lam_sum_chi2(lambda) = (sum(1 ./ (lambda[1] .+ beta.*eig_val))/size(eig_val,1)/size(q_BZ,1)-1/3)^2
    lam_lb = maximum(-beta.*eig_val)

    lower = lam_lb
    upper = Inf
    x0 = [lam_lb+10]
    result = optimize(lam_sum_chi2, lower, upper, x0)
    lambda = Optim.minimizer(result)[1]

    return lambda
end

function getCorr_iso(unitcell::UnitCell{D}, Jq_calc::Array{ComplexF64,3}, beta::Float64, lambda::Float64) where D
    natom = length(unitcell.basis)
    correl = zeros(axes(Jq_calc,3))
    Threads.@threads for ind_q in axes(Jq_calc,3)
        # correl[ind_q] = sum(inv(lambda.*Diagonal(ones(natom)) .+ beta.*Jq_calc[:,:,ind_q]));
        correl[ind_q] = real(sum(inv(lambda.*Matrix{Float64}(I,natom,natom) .+ beta.*Jq_calc[:,:,ind_q])));
    end
    return correl
end

function getChiT_iso(unitcell::UnitCell{D}, temperatures::Vector{Float64}) where D
    natom = length(unitcell.basis)
    beta_list = 1 ./ (temperatures ./ 11.6045)
    chiT::Vector{Float64} = zeros(axes(beta_list))
    for ind = 1:length(beta_list)
        beta = beta_list[ind]
        dist::Matrix{Float64} = getDist(unitcell)
        Jq_calc::Array{ComplexF64,3} = getFourier_iso(unitcell, dist, [0.0 0.0 0.0])
        lam_answer::Float64 = solveLambda_iso(unitcell, beta)
        chiT[ind] = getCorr_iso(unitcell, Jq_calc, beta, lam_answer)[1]/natom/3
    end
    return chiT
end


function getFourier_aniso(unitcell::UnitCell{D}, dist::Matrix{Float64}, q_crys::Matrix{Float64}) where D
    natom::Int64 = length(unitcell.basis)
    Jq::Array{ComplexF64, 3} = zeros(natom*3, natom*3, axes(q_crys,1))
    Threads.@threads for ind_q in axes(q_crys,1)
        q_lab = zeros(3)
        q_lab[1] = q_crys[ind_q,1]*unitcell.rl[1,1] + q_crys[ind_q,2]*unitcell.rl[2,1] + q_crys[ind_q,3]*unitcell.rl[3,1]
        q_lab[2] = q_crys[ind_q,1]*unitcell.rl[1,2] + q_crys[ind_q,2]*unitcell.rl[2,2] + q_crys[ind_q,3]*unitcell.rl[3,2]
        q_lab[3] = q_crys[ind_q,1]*unitcell.rl[1,3] + q_crys[ind_q,2]*unitcell.rl[2,3] + q_crys[ind_q,3]*unitcell.rl[3,3]

        nbond = length(unitcell.interactions)
        for ind_bond = 1:nbond
            bond = unitcell.interactions[ind_bond]
            phase::Float64 = q_lab[1]*dist[ind_bond, 1] + q_lab[2]*dist[ind_bond,2] + q_lab[3]*dist[ind_bond,3]
            Jterm::Matrix{ComplexF64} =  bond[4]*exp(-1im*phase) # no 1/2 factor, check logbook SCGA/A_further_check
            # Jq[3*(bond[1]-1)+1:3*(bond[1]-1)+3, 3*(bond[2]-1)+1:3*(bond[2]-1)+3, ind_q] += Jterm
            # Jq[3*(bond[2]-1)+1:3*(bond[2]-1)+3, 3*(bond[1]-1)+1:3*(bond[1]-1)+3, ind_q] += adjoint(Jterm)

            Jq[3*bond[1]-2, 3*bond[2]-2, ind_q] += Jterm[1]
            Jq[3*bond[1]-1, 3*bond[2]-2, ind_q] += Jterm[2]
            Jq[3*bond[1],   3*bond[2]-2, ind_q] += Jterm[3]
            Jq[3*bond[1]-2, 3*bond[2]-1, ind_q] += Jterm[4]
            Jq[3*bond[1]-1, 3*bond[2]-1, ind_q] += Jterm[5]
            Jq[3*bond[1],   3*bond[2]-1, ind_q] += Jterm[6]
            Jq[3*bond[1]-2, 3*bond[2],   ind_q] += Jterm[7]
            Jq[3*bond[1]-1, 3*bond[2],   ind_q] += Jterm[8]
            Jq[3*bond[1],   3*bond[2],   ind_q] += Jterm[9]

            Jq[3*bond[2]-2, 3*bond[1]-2, ind_q] += conj(Jterm[1])
            Jq[3*bond[2]-1, 3*bond[1]-2, ind_q] += conj(Jterm[4])
            Jq[3*bond[2],   3*bond[1]-2, ind_q] += conj(Jterm[7])
            Jq[3*bond[2]-2, 3*bond[1]-1, ind_q] += conj(Jterm[2])
            Jq[3*bond[2]-1, 3*bond[1]-1, ind_q] += conj(Jterm[5])
            Jq[3*bond[2],   3*bond[1]-1, ind_q] += conj(Jterm[8])
            Jq[3*bond[2]-2, 3*bond[1],   ind_q] += conj(Jterm[3])
            Jq[3*bond[2]-1, 3*bond[1],   ind_q] += conj(Jterm[6])
            Jq[3*bond[2],   3*bond[1],   ind_q] += conj(Jterm[9])
        end
        # interactionsOnsite
        for ind_atom = 1:natom
            Jq[3*ind_atom-2, 3*ind_atom-2, ind_q] += unitcell.interactionsOnsite[ind_atom][1,1]
            Jq[3*ind_atom-2, 3*ind_atom-1, ind_q] += unitcell.interactionsOnsite[ind_atom][1,2]
            Jq[3*ind_atom-2, 3*ind_atom,   ind_q] += unitcell.interactionsOnsite[ind_atom][1,3]

            Jq[3*ind_atom-1, 3*ind_atom-2, ind_q] += unitcell.interactionsOnsite[ind_atom][2,1]
            Jq[3*ind_atom-1, 3*ind_atom-1, ind_q] += unitcell.interactionsOnsite[ind_atom][2,2]
            Jq[3*ind_atom-1, 3*ind_atom,   ind_q] += unitcell.interactionsOnsite[ind_atom][2,3]

            Jq[3*ind_atom,   3*ind_atom-2, ind_q] += unitcell.interactionsOnsite[ind_atom][3,1]
            Jq[3*ind_atom,   3*ind_atom-1, ind_q] += unitcell.interactionsOnsite[ind_atom][3,2]
            Jq[3*ind_atom,   3*ind_atom,   ind_q] += unitcell.interactionsOnsite[ind_atom][3,3]
        end
    end
    return Jq
end

function solveLambda_aniso(unitcell::UnitCell{D}, beta::Float64) where D
    n_grid::Int64 = 8
    k0::Vector{Float64} = [(1-n_grid)/(2*n_grid):1/n_grid:(n_grid-1)/(2*n_grid);]
    qx::Array{Float64, 3} = [x_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    qy::Array{Float64, 3} = [y_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    qz::Array{Float64, 3} = [z_grid for x_grid in k0, y_grid in k0, z_grid in k0]
    q_BZ = [vec(qx) vec(qy) vec(qz)]

    dist = getDist(unitcell)
    Jq_BZ = getFourier_aniso(unitcell, dist, q_BZ)

    # eigensolution 
    natom = length(unitcell.basis)
    eig_val = zeros(natom*3, size(q_BZ,1))
    for ind_q in axes(q_BZ,1)
        eig_val[:,ind_q] = eigvals(Jq_BZ[:,:,ind_q])
    end

    lam_sum_chi2(lambda) = (sum(1 ./ (lambda[1] .+ beta.*eig_val))/size(eig_val,1)/size(q_BZ,1)-1/3)^2
    lam_lb = maximum(-beta.*eig_val)

    lower = lam_lb
    upper = Inf
    x0 = [lam_lb+10]
    result = optimize(lam_sum_chi2, lower, upper, x0)
    lambda = Optim.minimizer(result)[1]

    return lambda
end

function getCorr_aniso(unitcell::UnitCell{D}, q_crys::Matrix{Float64}, Jq_calc::Array{ComplexF64,3}, beta::Float64, lambda::Float64) where D
    correl::Vector{Float64} = zeros(axes(Jq_calc,3))
    Threads.@threads for ind_q in axes(q_crys,1)
        q_lab = zeros(3)
        q_lab[1] = q_crys[ind_q,1]*unitcell.rl[1,1] + q_crys[ind_q,2]*unitcell.rl[2,1] + q_crys[ind_q,3]*unitcell.rl[3,1]
        q_lab[2] = q_crys[ind_q,1]*unitcell.rl[1,2] + q_crys[ind_q,2]*unitcell.rl[2,2] + q_crys[ind_q,3]*unitcell.rl[3,2]
        q_lab[3] = q_crys[ind_q,1]*unitcell.rl[1,3] + q_crys[ind_q,2]*unitcell.rl[2,3] + q_crys[ind_q,3]*unitcell.rl[3,3]
        Q_norm::Vector{Float64} = q_lab/norm(q_lab)
        F = eigen(Jq_calc[:,:,ind_q])
        eig_val::Vector{Float64} = F.values
        U::Matrix{ComplexF64} = F.vectors

        for ind_eig in axes(U,2)
            Fsum = zeros(3)*im
            for ind_atom in axes(unitcell.basis,1)
                UQ::ComplexF64 = U[ind_atom*3-2, ind_eig]*Q_norm[1] + U[ind_atom*3-1, ind_eig]*Q_norm[2] + U[ind_atom*3, ind_eig]*Q_norm[3]
                Fsum[1] += U[ind_atom*3-2, ind_eig] - UQ*Q_norm[1]
                Fsum[2] += U[ind_atom*3-1, ind_eig] - UQ*Q_norm[2]
                Fsum[3] += U[ind_atom*3,   ind_eig] - UQ*Q_norm[3]
            end
            F2::Float64 = real(Fsum[1].*conj(Fsum[1]) + Fsum[2].*conj(Fsum[2]) + Fsum[3].*conj(Fsum[3]))
            denom::Float64 = lambda + beta*eig_val[ind_eig]
            correl[ind_q] += F2/denom
        end
    end
    return correl
end
