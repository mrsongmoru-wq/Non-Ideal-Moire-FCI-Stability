module WangContinuum

using LinearAlgebra

using ..TMBGZeroField

export AbstractWangModel,
       TBGWangParams,
       TBGWangModel,
       TMBGWangModel,
       build_tbg_wang,
       build_tmbg_wang,
       hamiltonian,
       central_eigensystem,
       gauge_shift_matrix,
       velocity_y,
       field_coupling_sign,
       reciprocal_vectors,
       direct_vectors,
       local_dimension,
       basis_dimension,
       magnetic_length,
       moire_area

abstract type AbstractWangModel end

const SIGMA0 = ComplexF64[1 0; 0 1]
const SIGMAX = ComplexF64[0 1; 1 0]
const SIGMAY = ComplexF64[0 -im; im 0]
const SIGMAZ = ComplexF64[1 0; 0 -1]

struct TBGWangParams
    theta_deg::Float64
    w0_ratio::Float64
    w1::Float64
    vf::Float64
end

TBGWangParams(; theta_deg=1.05, w0_ratio=0.7, w1=96.056, vf=2135.4) =
    TBGWangParams(theta_deg, w0_ratio, w1, vf)

struct TBGWangModel <: AbstractWangModel
    params::TBGWangParams
    lg::Int
    periodic_G::Bool
    g1::Vector{Float64}
    g2::Vector{Float64}
    a1::Vector{Float64}
    a2::Vector{Float64}
    Kt::Vector{Float64}
    Kb::Vector{Float64}
    Glist::Vector{Tuple{Int,Int}}
    Gdict::Dict{Tuple{Int,Int},Int}
    T0::Matrix{ComplexF64}
    T1::Matrix{ComplexF64}
    T2::Matrix{ComplexF64}
end

struct TMBGWangModel <: AbstractWangModel
    source_params::TMBGParams
    lg::Int
    periodic_G::Bool
    g1::Vector{Float64}
    g2::Vector{Float64}
    a1::Vector{Float64}
    a2::Vector{Float64}
    Glist::Vector{Tuple{Int,Int}}
    Gdict::Dict{Tuple{Int,Int},Int}
    kshift::Vector{Tuple{Float64,Float64}}
    vF::Float64
    t1::Float64
    v3::Float64
    v4::Float64
    hBN::Vector{Float64}
    Ts::Vector{Matrix{ComplexF64}}
end

function centered_square_grid(lg::Int)
    @assert isodd(lg) "Wang reciprocal grid requires odd lg"
    h = (lg - 1) ÷ 2
    glist = [(m, n) for n in -h:h for m in -h:h]
    return glist, Dict(g => i for (i, g) in enumerate(glist))
end

function real_lattice(g1::Vector{Float64}, g2::Vector{Float64})
    G = hcat(g1, g2)
    A = 2pi .* inv(G)'
    return Vector(A[:, 1]), Vector(A[:, 2])
end

function build_tbg_wang(params::TBGWangParams=TBGWangParams(); lg::Int=9, periodic_G::Bool=true)
    theta = deg2rad(params.theta_deg)
    kb = 8pi / 3 * sin(theta / 2)
    gmag = sqrt(3) * kb
    g1 = [gmag, 0.0]
    g2 = gmag .* [cos(2pi / 3), sin(2pi / 3)]
    a1, a2 = real_lattice(g1, g2)
    Kt = kb .* [cos(5pi / 6), sin(5pi / 6)]
    Kb = kb .* [cos(-5pi / 6), sin(-5pi / 6)]
    Glist, Gdict = centered_square_grid(lg)

    w0 = params.w0_ratio * params.w1
    omega = exp(2pi * im / 3)
    T0 = ComplexF64[w0 params.w1; params.w1 w0]
    T1 = ComplexF64[w0 params.w1 * conj(omega); params.w1 * omega w0]
    T2 = ComplexF64[w0 params.w1 * omega; params.w1 * conj(omega) w0]
    return TBGWangModel(params, lg, periodic_G, g1, g2, a1, a2, Kt, Kb,
                        Glist, Gdict, T0, T1, T2)
end

function build_tmbg_wang(params::TMBGParams; lg::Int=9, periodic_G::Bool=true)
    source = build_model(params)
    # Unimodular reciprocal-basis change used by Wang and Vafek:
    # g1W = g1 - g2, g2W = g2. The corresponding direct basis is a1W=a1,
    # a2W=a1+a2 and has the 60 degree convention used in the MTG derivation.
    g1 = source.g1 .- source.g2
    g2 = copy(source.g2)
    a1, a2 = real_lattice(g1, g2)
    Glist, Gdict = centered_square_grid(lg)
    kshift = [(1 / 3, 2 / 3), (0.0, 0.0), (0.0, 0.0)]
    return TMBGWangModel(params, lg, periodic_G, g1, g2, a1, a2, Glist, Gdict,
                         kshift, source.vF, source.t1, source.v3, source.v4,
                         source.hBN, source.Ts)
end

local_dimension(::TBGWangModel) = 4
local_dimension(::TMBGWangModel) = 6
basis_dimension(model::AbstractWangModel) = local_dimension(model) * length(model.Glist)
reciprocal_vectors(model::AbstractWangModel) = (model.g1, model.g2)
direct_vectors(model::AbstractWangModel) = (model.a1, model.a2)
moire_area(model::AbstractWangModel) = abs(det(hcat(model.a1, model.a2)))
magnetic_length(model::AbstractWangModel, p::Int, q::Int) =
    sqrt(moire_area(model) * q / (2pi * abs(p)))

function wrapped_key(model::AbstractWangModel, m::Int, n::Int; periodic::Bool=model.periodic_G)
    if !periodic
        return haskey(model.Gdict, (m, n)) ? (m, n) : nothing
    end
    h = (model.lg - 1) ÷ 2
    return (mod(m + h, model.lg) - h, mod(n + h, model.lg) - h)
end

function shifted_index(model::AbstractWangModel, g::Tuple{Int,Int}, shift::Tuple{Int,Int};
                       periodic::Bool=model.periodic_G)
    key = wrapped_key(model, g[1] + shift[1], g[2] + shift[2]; periodic=periodic)
    return key === nothing ? nothing : model.Gdict[key]
end

function gauge_shift_matrix(model::AbstractWangModel, dn1::Int, dn2::Int; periodic::Bool=true)
    nlocal = local_dimension(model)
    nG = length(model.Glist)
    V = zeros(ComplexF64, nlocal * nG, nlocal * nG)
    for (i, g) in enumerate(model.Glist)
        j = shifted_index(model, g, (dn1, dn2); periodic=periodic)
        j === nothing && continue
        rr = (nlocal * (i - 1) + 1):(nlocal * i)
        cc = (nlocal * (j - 1) + 1):(nlocal * j)
        V[rr, cc] .= I(nlocal)
    end
    return V
end

@inline function dirac_block(k::Vector{Float64}, vf::Float64)
    return vf .* ComplexF64[0 k[1] - im * k[2]; k[1] + im * k[2] 0]
end

function add_tunneling!(H::Matrix{ComplexF64}, model::TBGWangModel,
                        i::Int, j::Union{Nothing,Int}, T::Matrix{ComplexF64})
    j === nothing && return
    ri = (4 * (i - 1) + 1):(4 * (i - 1) + 2)
    rj = (4 * (j - 1) + 3):(4 * j)
    H[ri, rj] .+= T
    H[rj, ri] .+= T'
end

function hamiltonian(model::TBGWangModel, k::Vector{Float64})
    H = zeros(ComplexF64, basis_dimension(model), basis_dimension(model))
    for (i, g) in enumerate(model.Glist)
        G = g[1] .* model.g1 .+ g[2] .* model.g2
        r1 = (4 * (i - 1) + 1):(4 * (i - 1) + 2)
        r2 = (4 * (i - 1) + 3):(4 * i)
        H[r1, r1] .+= dirac_block(k .- model.Kb .+ G, model.params.vf)
        H[r2, r2] .+= dirac_block(k .- model.Kt .+ G, model.params.vf)

        add_tunneling!(H, model, i, shifted_index(model, g, (0, 0)), model.T0)
        add_tunneling!(H, model, i, shifted_index(model, g, (0, 1)), model.T2)
        add_tunneling!(H, model, i, shifted_index(model, g, (1, 1)), model.T1)
    end
    return H
end

function tmbg_local_block(model::TMBGWangModel, k::Vector{Float64}, g::Tuple{Int,Int})
    H = zeros(ComplexF64, 6, 6)
    Udiag = [-model.source_params.U / 2, 0.0, model.source_params.U / 2]
    H .+= kron(Diagonal(model.hBN), SIGMAZ)
    H .+= kron(Diagonal(Udiag), SIGMA0)
    G = g[1] .* model.g1 .+ g[2] .* model.g2

    kvals = Vector{Float64}[]
    for layer in 1:3
        sx, sy = model.kshift[layer]
        push!(kvals, k .+ G .+ sx .* model.g1 .+ sy .* model.g2)
        d = model.source_params.valley * kvals[end][1] - im * kvals[end][2]
        rr = (2 * layer - 1):(2 * layer)
        H[rr, rr] .+= -model.vF .* ComplexF64[0 d; conj(d) 0]
    end

    kval = kvals[2]
    kp = model.source_params.valley * kval[1] + im * kval[2]
    km = model.source_params.valley * kval[1] - im * kval[2]
    HAB = ComplexF64[model.v3 * kp model.t1;
                     model.v4 * km model.v3 * kp]
    H[5:6, 3:4] .+= HAB
    H[3:4, 5:6] .+= HAB'
    return H
end

function add_tmbg_tunneling!(H::Matrix{ComplexF64}, model::TMBGWangModel,
                             i::Int, j::Union{Nothing,Int}, T::Matrix{ComplexF64})
    j === nothing && return
    ri = (6 * (i - 1) + 1):(6 * (i - 1) + 2)
    rj = (6 * (j - 1) + 3):(6 * (j - 1) + 4)
    H[ri, rj] .+= T
    H[rj, ri] .+= T'
end

function hamiltonian(model::TMBGWangModel, k::Vector{Float64})
    H = zeros(ComplexF64, basis_dimension(model), basis_dimension(model))
    for (i, g) in enumerate(model.Glist)
        rr = (6 * (i - 1) + 1):(6 * i)
        H[rr, rr] .+= tmbg_local_block(model, k, g)
        add_tmbg_tunneling!(H, model, i, shifted_index(model, g, (0, 0)), model.Ts[1])
        add_tmbg_tunneling!(H, model, i, shifted_index(model, g, (1, 1)), model.Ts[2])
        add_tmbg_tunneling!(H, model, i, shifted_index(model, g, (0, 1)), model.Ts[3])
    end
    return H
end

function velocity_y(model::TBGWangModel)
    return model.params.vf .* kron(ComplexF64[1 0; 0 1], SIGMAY)
end

function velocity_y(model::TMBGWangModel)
    tau2 = zeros(ComplexF64, 3, 3)
    tau2[3, 2] = 1
    D = ComplexF64[im * model.v3 0; -im * model.v4 im * model.v3]
    return -model.vF .* kron(Matrix{ComplexF64}(I, 3, 3), SIGMAY) .+
           kron(tau2, D) .+ kron(tau2', D')
end

# With the reciprocal-index and MTG-phase orientation used below, positive
# signB implements k_y -> k_y - x/l_B^2. The field direction must be reversed
# together with every MTG phase, not by changing this prefactor alone.
field_coupling_sign(::TBGWangModel, signB::Int) = -signB
field_coupling_sign(::TMBGWangModel, signB::Int) = -signB

function central_eigensystem(model::AbstractWangModel, k::Vector{Float64}, nactive::Int)
    @assert iseven(nactive)
    H = hamiltonian(model, k)
    nb = size(H, 1)
    start = nb ÷ 2 - nactive ÷ 2 + 1
    ids = start:(start + nactive - 1)

    # The tMBG central bands are not separated by the sign of their energy:
    # one member of the target pair crosses zero across the Brillouin zone.
    # LAPACK's index-range driver preserves the global band indices while
    # computing only the requested eigenpairs.
    F = eigen(Hermitian(H), ids)
    return F.values, F.vectors
end

end
