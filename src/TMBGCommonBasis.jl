module TMBGCommonBasis

using LinearAlgebra

using ..TMBGMagneticHW
using ..TMBGMagneticHW.TMBGHybridWannier.WangContinuum

export CommonBasisProjection,
       reciprocal_shift_map,
       lifted_row_amplitudes!,
       inverse_sqrt_metric

struct CommonBasisProjection
    selected_bands::Vector{Int}
    band_weights::Array{Float64,3}
    target_weights::Matrix{Float64}
    active_weights::Matrix{Float64}
    common_weights::Matrix{Float64}
    source_weights::Vector{Float64}
    target_fraction_within_active::Float64
    raw_metric_errors::Vector{Float64}
    raw_metric_eigenvalue_bounds::Matrix{Float64}
end

@inline trial_index(alpha::Int, ir::Int, in_index::Int, q::Int, nactive::Int) =
    alpha + nactive * (ir - 1) + nactive * q * (in_index - 1)

function reciprocal_shift_map(
    model::AbstractWangModel,
    dg1::Int,
    dg2::Int;
    periodic::Bool=model.periodic_G,
)
    half_width = div(model.lg - 1, 2)
    result = zeros(Int, length(model.Glist))
    for (index, reciprocal) in enumerate(model.Glist)
        key = (reciprocal[1] + dg1, reciprocal[2] + dg2)
        if periodic
            key = (
                mod(key[1] + half_width, model.lg) - half_width,
                mod(key[2] + half_width, model.lg) - half_width,
            )
        end
        result[index] = get(model.Gdict, key, 0)
    end
    return result
end

function add_shifted!(
    destination::Matrix{ComplexF64},
    source::Matrix{ComplexF64},
    shift_map::Vector{Int},
    nlocal::Int,
    phase::ComplexF64,
)
    for reciprocal_index in eachindex(shift_map)
        shifted_index = shift_map[reciprocal_index]
        shifted_index == 0 && continue
        rows = (nlocal * (reciprocal_index - 1) + 1):(nlocal * reciprocal_index)
        columns = (nlocal * (shifted_index - 1) + 1):(nlocal * shifted_index)
        @views destination[rows, :] .+= phase .* source[columns, :]
    end
    return destination
end

function lifted_row_amplitudes!(
    amplitudes::Matrix{ComplexF64},
    result::MagneticHWModel,
    coefficients::Matrix{ComplexF64},
    external_i1::Int,
    i2::Int,
    i1_row::Int,
    ir_row::Int,
    translations::UnitRange{Int},
    shift_maps::Dict{Tuple{Int,Int},Vector{Int}},
)
    fill!(amplitudes, 0)
    model, hw, flux = result.model, result.hw, result.flux
    nactive = size(hw.WLS, 2)
    nstates = size(coefficients, 2)
    q, l1, l2, nq = flux.q, result.l1, result.l2, result.nq
    p_signed = flux.signB * flux.p
    nlocal = local_dimension(model)
    source = zeros(ComplexF64, nactive, nstates)
    bloch = zeros(ComplexF64, basis_dimension(model), nstates)
    external_k1 = result.k1[external_i1]

    for translation in translations
        i1_unwrapped = i1_row - translation * p_signed * nq
        ir_unwrapped = ir_row - translation * p_signed
        i1_column = mod1(i1_unwrapped, l1)
        ir_column = mod1(ir_unwrapped, q)
        dg1 = div(i1_unwrapped - i1_column, l1)
        dg2 = div(ir_unwrapped - ir_column, q)
        shift_map = shift_maps[(dg1, dg2)]
        i2_column = (ir_column - 1) * l2 + i2
        reduced_k1 = (i1_column - 1) / l1

        fill!(source, 0)
        for (center_index, center) in enumerate(result.nvec)
            rows = [
                trial_index(alpha, ir_column, center_index, q, nactive)
                for alpha in 1:nactive
            ]
            @views source .+=
                exp(-im * 2pi * reduced_k1 * center) .* coefficients[rows, :]
        end
        mul!(bloch, @view(hw.WLS[:, :, i1_column, i2_column]), source)
        phase = exp(im * 2pi * external_k1 * translation) *
                exp(-im * pi * translation * (translation - 1) * p_signed / (2q)) *
                exp(-im * 2pi * reduced_k1 * translation) / l1
        add_shifted!(amplitudes, bloch, shift_map, nlocal, phase)
    end
    return amplitudes
end

function inverse_sqrt_metric(
    metric::Hermitian{ComplexF64,Matrix{ComplexF64}};
    tolerance::Float64=1e-12,
)
    factorization = eigen(metric)
    minimum(factorization.values) > tolerance ||
        error("Lifted common-basis subspace is rank deficient")
    return factorization.vectors *
           Diagonal(1 ./ sqrt.(factorization.values)) *
           factorization.vectors'
end

end
