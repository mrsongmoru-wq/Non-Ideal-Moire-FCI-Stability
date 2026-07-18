module TMBGZeroField

using LinearAlgebra

export TMBGParams, TMBGModel, build_model, hamiltonian, band_path, solve_band_path

struct TMBGParams
    theta_deg::Float64
    N::Int
    kappa::Float64
    valley::Int
    hBN_potential::Float64
    hBNlayer::NTuple{3,Int}
    U::Float64
end

mutable struct TMBGModel
    params::TMBGParams
    vF::Float64
    t1::Float64
    g3::Float64
    g4::Float64
    v3::Float64
    v4::Float64
    w::Float64
    ktheta::Float64
    hBN::Vector{Float64}
    g1::Vector{Float64}
    g2::Vector{Float64}
    kshift::Vector{Tuple{Float64,Float64}}
    Glist::Vector{Tuple{Int,Int}}
    Gdict::Dict{Tuple{Int,Int},Int}
    sizeH_h::Int
    Ts::Vector{Matrix{ComplexF64}}
    t_s::Vector{Matrix{ComplexF64}}
end

const A_REF = 2.46
const W1_DEFAULT = 110.0

const SIGMA0 = ComplexF64[1 0; 0 1]
const SIGMAX = ComplexF64[0 1; 1 0]
const SIGMAY = ComplexF64[0 -im; im 0]
const SIGMAZ = ComplexF64[1 0; 0 -1]
const SIGMAP = (SIGMAX + im * SIGMAY) / 2

function tau1()
    t = zeros(ComplexF64, 3, 3)
    t[1, 2] = 1
    return t
end

function tau2()
    t = zeros(ComplexF64, 3, 3)
    t[3, 2] = 1
    return t
end

function build_model(params::TMBGParams)
    vF = sqrt(3) / 2 * A_REF * 2610.0
    t1 = 361.0
    g3 = 140.0
    g4 = 283.0
    v3 = sqrt(3) / 2 * A_REF * g3
    v4 = sqrt(3) / 2 * A_REF * g4
    w = W1_DEFAULT
    ktheta = 8 * pi / 3 / A_REF * sin(params.theta_deg * pi / 360)
    hBN = params.hBN_potential / 2 .* collect(Float64, params.hBNlayer)
    g1 = sqrt(3) * ktheta .* [0.5, sqrt(3) / 2]
    g2 = sqrt(3) * ktheta .* [-0.5, sqrt(3) / 2]
    kshift = [(1 / 3, 1 / 3), (0.0, 0.0), (0.0, 0.0)]

    N = params.N
    Glist = Tuple{Int,Int}[]
    for m in -N:N, n in -N:N
        if abs(m + n) <= N
            push!(Glist, (m, n))
        end
    end
    Gdict = Dict(g => i for (i, g) in enumerate(Glist))
    sizeH_h = length(Glist) * 3

    Ts = Matrix{ComplexF64}[]
    for s in 0:2
        push!(Ts, w .* (params.kappa .* SIGMA0 .+
                        cos(params.valley * 2pi * s / 3) .* SIGMAX .+
                        sin(params.valley * 2pi * s / 3) .* SIGMAY))
    end
    t_s = [kron(tau1(), T) for T in Ts]

    return TMBGModel(params, vF, t1, g3, g4, v3, v4, w, ktheta, hBN,
                     g1, g2, kshift, Glist, Gdict, sizeH_h, Ts, t_s)
end

function ab_layer_cp(model::TMBGModel, k::Vector{Float64})
    valley = model.params.valley
    kplus = valley * k[1] + im * k[2]
    kminus = valley * k[1] - im * k[2]
    HAB = ComplexF64[model.v3 * kplus model.t1;
                     model.v4 * kminus model.v3 * kplus]
    return kron(tau2(), HAB)
end

function hamiltonian(model::TMBGModel, k_reduce::Vector{Float64})
    nG = length(model.Glist)
    Hsize = 6 * nG
    H0 = zeros(ComplexF64, Hsize, Hsize)
    HT = zeros(ComplexF64, Hsize, Hsize)
    Udiag = [-model.params.U / 2, 0.0, model.params.U / 2]
    hBNmat = kron(Diagonal(model.hBN), SIGMAZ) / 2
    Umat = kron(Diagonal(Udiag), SIGMA0) / 2

    for (ig, (n1, n2)) in enumerate(model.Glist)
        r = (6 * (ig - 1) + 1):(6 * ig)
        H0[r, r] .+= hBNmat
        H0[r, r] .+= Umat

        kvals = Vector{Float64}[]
        for L in 1:3
            sx, sy = model.kshift[L]
            push!(kvals, k_reduce .+ (n1 + sx) .* model.g1 .+ (n2 + sy) .* model.g2)
        end
        diag_vals = ComplexF64[
            model.params.valley * kvals[L][1] - im * kvals[L][2] for L in 1:3
        ]
        H0[r, r] .+= -model.vF .* kron(Diagonal(diag_vals), SIGMAP)
        H0[r, r] .+= ab_layer_cp(model, kvals[2])

        HT[r, r] .+= model.t_s[1]
        for (shift, ts_index) in (((1, 0), 2), ((0, 1), 3))
            key = (n1 + shift[1], n2 + shift[2])
            if haskey(model.Gdict, key)
                jg = model.Gdict[key]
                c = (6 * (jg - 1) + 1):(6 * jg)
                HT[r, c] .+= model.t_s[ts_index]
            end
        end
    end
    H = H0 + HT
    return H + H'
end

function path_points(num_total_points::Int; mode::Int=0)
    Gamma = [2 / 3, -1 / 3]
    Gamma2 = [-1 / 3, -1 / 3]
    K1 = [1 / 3, 1 / 3]
    K2 = [0.0, 0.0]
    M = (K1 .+ K2) ./ 2
    if mode == 0
        labels = ["K", "Gamma", "M", "Kprime"]
        kpath = [K1, Gamma, M, K2]
    elseif mode == 1
        labels = ["Kprime", "K", "Gamma2", "Gamma", "Kprime"]
        kpath = [K2, K1, Gamma2, Gamma, K2]
    else
        error("Unsupported path mode $mode")
    end

    return labels, kpath
end

function band_path(model::TMBGModel, num_total_points::Int; mode::Int=0)
    labels, kpath = path_points(num_total_points; mode=mode)
    distances = Float64[]
    for i in 1:(length(kpath)-1)
        dk = (kpath[i+1][1] - kpath[i][1]) .* model.g1 .+
             (kpath[i+1][2] - kpath[i][2]) .* model.g2
        push!(distances, norm(dk))
    end
    total = sum(distances)
    npts = [max(2, round(Int, d / total * num_total_points)) for d in distances]
    diff = num_total_points - sum(npts)
    while diff != 0
        for i in eachindex(npts)
            diff == 0 && break
            if diff > 0
                npts[i] += 1
                diff -= 1
            elseif npts[i] > 2
                npts[i] -= 1
                diff += 1
            end
        end
    end

    points = Vector{Vector{Float64}}()
    segment_lengths = [0]
    total_points = 0
    for i in 1:(length(kpath)-1)
        nseg = npts[i]
        endpoint = (i == length(kpath)-1)
        count = endpoint ? nseg - 1 : nseg
        for j in 0:(count-1)
            t = endpoint && count > 1 ? j / (count - 1) : j / nseg
            push!(points, (1 - t) .* kpath[i] .+ t .* kpath[i+1])
        end
        total_points += count
        push!(segment_lengths, total_points)
    end
    return labels, segment_lengths, points
end

function solve_band_path(model::TMBGModel, num_total_points::Int; mode::Int=0)
    labels, segment_lengths, points = band_path(model, num_total_points; mode=mode)
    nb = 2 * model.sizeH_h
    energies = zeros(Float64, length(points), nb)
    for (i, kred) in enumerate(points)
        kcart = kred[1] .* model.g1 .+ kred[2] .* model.g2
        vals = eigvals(Hermitian(hamiltonian(model, kcart)))
        energies[i, :] .= vals
    end
    return labels, segment_lengths, points, energies
end

end
