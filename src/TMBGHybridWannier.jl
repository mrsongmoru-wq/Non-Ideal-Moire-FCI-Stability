module TMBGHybridWannier

using LinearAlgebra

include("TMBGZeroField.jl")
using .TMBGZeroField
include("WangContinuum.jl")
using .WangContinuum

export HybridWannierTMBG,
       build_hybrid_wannier,
       energy_eigenbasis,
       hybrid_to_energy_transform,
       flat_band_grid,
       polar_unitary,
       TMBGZeroField,
       WangContinuum

mutable struct HybridWannierTMBG
    model::AbstractWangModel
    n1::Int
    n2::Int
    band_indices::Vector{Int}
    eigvals::Array{Float64,3}
    eigvecs::Array{ComplexF64,4}
    bloch_ham::Array{ComplexF64,4}
    links::Array{ComplexF64,4}
    partial::Array{ComplexF64,4}
    wilson::Array{ComplexF64,3}
    wilson_vals::Array{ComplexF64,2}
    wilson_vecs::Array{ComplexF64,3}
    coeffs::Array{ComplexF64,4}
    WLS::Array{ComplexF64,4}
    Hwk::Array{ComplexF64,4}
    chern::Int
    chern_raw::Float64
end

function energy_eigenbasis(hw::HybridWannierTMBG, i1::Int, i2::Int)
    frame_hamiltonian = Hermitian(Matrix(@view hw.bloch_ham[:, :, i1, i2]))
    factorization = eigen(frame_hamiltonian)
    frame = @view hw.eigvecs[:, :, i1, i2]
    vectors = Matrix(frame * factorization.vectors)
    return factorization.values, vectors
end

function hybrid_to_energy_transform(hw::HybridWannierTMBG, i1::Int, i2::Int)
    _, energy_vectors = energy_eigenbasis(hw, i1, i2)
    wannier_frame = @view hw.WLS[:, :, i1, i2]
    return Matrix(energy_vectors' * wannier_frame)
end

function polar_unitary(M::AbstractMatrix{ComplexF64})
    F = svd(M)
    return F.U * F.Vt
end

function nonabelian_chern(model::AbstractWangModel, vecs::Array{ComplexF64,4})
    _, _, n1, n2 = size(vecs)
    link1 = zeros(ComplexF64, n1, n2)
    link2 = zeros(ComplexF64, n1, n2)
    sewing1 = gauge_shift_matrix(model, 1, 0; periodic=true)
    sewing2 = gauge_shift_matrix(model, 0, 1; periodic=true)

    for i2 in 1:n2, i1 in 1:n1
        u = @view vecs[:, :, i1, i2]
        v1 = i1 < n1 ? Matrix(@view vecs[:, :, i1 + 1, i2]) :
             sewing1 * Matrix(@view vecs[:, :, 1, i2])
        v2 = i2 < n2 ? Matrix(@view vecs[:, :, i1, i2 + 1]) :
             sewing2 * Matrix(@view vecs[:, :, i1, 1])
        link1[i1, i2] = det(polar_unitary(Matrix(u' * v1)))
        link2[i1, i2] = det(polar_unitary(Matrix(u' * v2)))
    end

    flux = 0.0
    for i2 in 1:n2, i1 in 1:n1
        ip = mod1(i1 + 1, n1)
        jp = mod1(i2 + 1, n2)
        loop = link1[i1, i2] * link2[ip, i2] *
               conj(link1[i1, jp]) * conj(link2[i1, i2])
        flux += angle(loop)
    end
    raw = flux / (2pi)
    integer = round(Int, raw)
    abs(raw - integer) < 1e-6 || error("Active-subspace Chern number is not converged: $raw")
    return integer, raw
end

function central_band_indices(model::AbstractWangModel, nactive::Int)
    nb = basis_dimension(model)
    start = nb ÷ 2 - nactive ÷ 2 + 1
    return collect(start:(start + nactive - 1))
end

function flat_band_grid(model::AbstractWangModel, n1::Int, n2::Int=n1; nactive::Int=2)
    nbasis = basis_dimension(model)
    vals = zeros(Float64, nactive, n1, n2)
    vecs = zeros(ComplexF64, nbasis, nactive, n1, n2)
    hactive = zeros(ComplexF64, nactive, nactive, n1, n2)
    for i2 in 1:n2, i1 in 1:n1
        kred1 = (i1 - 1) / n1
        kred2 = (i2 - 1) / n2
        kcart = kred1 .* model.g1 .+ kred2 .* model.g2
        evals, evecs = central_eigensystem(model, kcart, nactive)
        vals[:, i1, i2] .= evals
        vecs[:, :, i1, i2] .= evecs
        hactive[:, :, i1, i2] .= Diagonal(evals)
    end
    return vals, vecs, hactive, central_band_indices(model, nactive)
end

function flat_band_grid(model::TMBGModel, n1::Int, n2::Int=n1; nactive::Int=2, lg::Int=2model.params.N + 1)
    return flat_band_grid(build_tmbg_wang(model.params; lg=lg), n1, n2; nactive=nactive)
end

function fix_k2_seed_gauge!(vecs::Array{ComplexF64,4},
                            hactive::Array{ComplexF64,4}, n2::Int)
    for i2 in 1:(n2 - 1)
        u = @view vecs[:, :, 1, i2]
        v = @view vecs[:, :, 1, i2 + 1]
        R = polar_unitary(Matrix(u' * v))'
        v[:, :] .= v * R
        h = Matrix(@view hactive[:, :, 1, i2 + 1])
        hactive[:, :, 1, i2 + 1] .= R' * h * R
    end
end

function best_permutation(weights::Matrix{Float64})
    n = size(weights, 1)
    best = collect(1:n)
    current = zeros(Int, n)
    used = falses(n)
    best_score = Ref(-Inf)

    function search(row::Int, score::Float64)
        if row > n
            if score > best_score[]
                best_score[] = score
                best .= current
            end
            return
        end
        for col in 1:n
            used[col] && continue
            used[col] = true
            current[row] = col
            search(row + 1, score + weights[row, col])
            used[col] = false
        end
    end
    search(1, 0.0)
    return best
end

function ordered_wilson_eigensystem(W::Matrix{ComplexF64}, previous)
    F = eigen(W)
    if previous === nothing
        order = sortperm(mod.(angle.(F.values), 2pi))
        return F.values[order], F.vectors[:, order]
    end

    order = best_permutation(abs.(previous' * F.vectors))
    vals = F.values[order]
    vecs = F.vectors[:, order]
    for a in axes(vecs, 2)
        z = dot(previous[:, a], vecs[:, a])
        abs(z) > 1e-14 && (vecs[:, a] .*= exp(-im * angle(z)))
    end
    return vals, vecs
end

function build_hybrid_wannier(model::AbstractWangModel, n1::Int, n2::Int=n1;
                              nactive::Int=2)
    vals, vecs, hactive, ids = flat_band_grid(model, n1, n2; nactive=nactive)
    chern, chern_raw = nonabelian_chern(model, vecs)
    fix_k2_seed_gauge!(vecs, hactive, n2)
    nbasis = size(vecs, 1)

    links = zeros(ComplexF64, nactive, nactive, n1, n2)
    partial = zeros(ComplexF64, nactive, nactive, n1, n2)
    wilson = zeros(ComplexF64, nactive, nactive, n2)
    wilson_vals = zeros(ComplexF64, nactive, n2)
    wilson_vecs = zeros(ComplexF64, nactive, nactive, n2)
    coeffs = zeros(ComplexF64, nactive, nactive, n1, n2)
    WLS = zeros(ComplexF64, nbasis, nactive, n1, n2)
    Hwk = zeros(ComplexF64, nactive, nactive, n1, n2)
    sewing = gauge_shift_matrix(model, 1, 0; periodic=true)
    previous = nothing

    for i2 in 1:n2
        for i1 in 1:n1
            u = @view vecs[:, :, i1, i2]
            v = i1 < n1 ? (@view vecs[:, :, i1 + 1, i2]) : sewing * (@view vecs[:, :, 1, i2])
            links[:, :, i1, i2] .= polar_unitary(Matrix(u' * v))
        end

        partial[:, :, 1, i2] .= I(nactive)
        for i1 in 2:n1
            partial[:, :, i1, i2] .= partial[:, :, i1 - 1, i2] * links[:, :, i1 - 1, i2]
        end
        W = partial[:, :, n1, i2] * links[:, :, n1, i2]
        wilson[:, :, i2] .= W
        wvals, wvecs = ordered_wilson_eigensystem(Matrix(W), previous)
        wilson_vals[:, i2] .= wvals
        wilson_vecs[:, :, i2] .= wvecs
        previous = copy(wvecs)

        for i1 in 1:n1
            phase = Diagonal(wvals .^ ((i1 - 1) / n1))
            C = partial[:, :, i1, i2]' * wvecs * phase
            coeffs[:, :, i1, i2] .= C
            u = @view vecs[:, :, i1, i2]
            WLS[:, :, i1, i2] .= u * C
            h = @view hactive[:, :, i1, i2]
            Hwk[:, :, i1, i2] .= C' * h * C
        end
    end

    return HybridWannierTMBG(model, n1, n2, ids, vals, vecs, hactive, links,
                             partial, wilson, wilson_vals, wilson_vecs, coeffs,
                             WLS, Hwk, chern, chern_raw)
end

function build_hybrid_wannier(model::TMBGModel, n1::Int, n2::Int=n1;
                              nactive::Int=2, lg::Int=2model.params.N + 1)
    wang_model = build_tmbg_wang(model.params; lg=lg)
    return build_hybrid_wannier(wang_model, n1, n2; nactive=nactive)
end

end
