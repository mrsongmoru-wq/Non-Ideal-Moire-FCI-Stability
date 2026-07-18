module TMBGMagneticGeometry

using FFTW
using LinearAlgebra

using ..TMBGMagneticHW
using ..TMBGMagneticHW.TMBGHybridWannier.WangContinuum
using ..TMBGCommonBasis

export MagneticGeometryResult,
       lift_magnetic_real_space_states,
       cartesian_metric_from_projector_distances,
       compute_magnetic_geometry

struct MagneticGeometryResult
    band::Int
    berry_flux::Matrix{Float64}
    berry_curvature::Matrix{Float64}
    metric_trace::Matrix{Float64}
    trace_violation::Matrix{Float64}
    metric_xx::Matrix{Float64}
    metric_yy::Matrix{Float64}
    metric_xy::Matrix{Float64}
    link_x_abs::Matrix{Float64}
    link_y_abs::Matrix{Float64}
    chern::Float64
    eta::Float64
    integrated_trace::Float64
    integrated_abs_curvature::Float64
    plaquette_area::Float64
end

function magnetic_shift_maps(
    result::MagneticHWModel,
    translations::UnitRange{Int};
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
)
    l1, q, nq = result.l1, result.flux.q, result.nq
    p_signed = result.flux.signB * result.flux.p
    shifts = Set{Tuple{Int,Int}}()
    for ir_row in 1:q, i1_row in 1:l1, translation in translations
        i1_unwrapped = i1_row - translation * p_signed * nq
        ir_unwrapped = ir_row - translation * p_signed
        push!(shifts, (
            div(i1_unwrapped - mod1(i1_unwrapped, l1), l1),
            div(ir_unwrapped - mod1(ir_unwrapped, q), q),
        ))
    end
    return Dict(
        shift => TMBGCommonBasis.reciprocal_shift_map(
            result.model,
            shift...;
            periodic=reciprocal_shift_periodic,
        )
        for shift in shifts
    )
end

function lift_magnetic_momentum_state(
    result::MagneticHWModel,
    band::Int,
    external_i1::Int,
    i2::Int,
    translations::UnitRange{Int},
    shift_maps::Dict{Tuple{Int,Int},Vector{Int}},
)
    model = result.model
    l1, l2, q = result.l1, result.l2, result.flux.q
    nlocal = local_dimension(model)
    fine_dimension1 = l1 * model.lg
    fine_dimension2 = l2 * q * model.lg
    ik = external_i1 + (i2 - 1) * l1
    coefficients = result.uort[:, :, ik] *
                   @view(result.eigvecs[:, band:band, ik])
    plane_wave_amplitudes = zeros(ComplexF64, basis_dimension(model), 1)
    fine_momentum = zeros(ComplexF32, nlocal, fine_dimension1, fine_dimension2)

    for ir_row in 1:q, i1_row in 1:l1
        TMBGCommonBasis.lifted_row_amplitudes!(
            plane_wave_amplitudes,
            result,
            coefficients,
            external_i1,
            i2,
            i1_row,
            ir_row,
            translations,
            shift_maps,
        )
        i2_full = (ir_row - 1) * l2 + i2
        for (g_index, g) in enumerate(model.Glist)
            fine_i1 = mod((i1_row - 1) + l1 * g[1], fine_dimension1) + 1
            fine_i2 = mod((i2_full - 1) + l2 * q * g[2], fine_dimension2) + 1
            rows = (nlocal * (g_index - 1) + 1):(nlocal * g_index)
            @views fine_momentum[:, fine_i1, fine_i2] .=
                ComplexF32.(plane_wave_amplitudes[rows, 1])
        end
    end
    normalization = sqrt(sum(abs2, fine_momentum))
    normalization > 1e-12 || error("Lifted magnetic state has zero norm")
    fine_momentum ./= normalization
    return fine_momentum
end

function lift_magnetic_momentum_row(
    result::MagneticHWModel,
    band::Int,
    i2::Int,
    translations::UnitRange{Int},
    shift_maps::Dict{Tuple{Int,Int},Vector{Int}},
)
    model = result.model
    fine_dimension1 = result.l1 * model.lg
    fine_dimension2 = result.l2 * result.flux.q * model.lg
    nlocal = local_dimension(model)
    row = Array{ComplexF32}(
        undef,
        nlocal,
        fine_dimension1,
        fine_dimension2,
        result.l1,
    )
    for i1 in 1:result.l1
        @views row[:, :, :, i1] .= lift_magnetic_momentum_state(
            result,
            band,
            i1,
            i2,
            translations,
            shift_maps,
        )
    end
    return row
end

function transported_overlap_momentum(
    bra::AbstractArray{<:Complex,3},
    ket::AbstractArray{<:Complex,3},
    shift1::Int,
    shift2::Int,
)
    dimension1, dimension2 = size(bra, 2), size(bra, 3)
    value = 0.0 + 0.0im
    # Multiplication by exp[-i delta.k r] moves the ket coefficient at
    # momentum index (j1+shift1,j2+shift2) to the bra index (j1,j2).
    for j2 in 1:dimension2
        source_j2 = mod1(j2 + shift2, dimension2)
        for j1 in 1:dimension1
            source_j1 = mod1(j1 + shift1, dimension1)
            for local_index in axes(bra, 1)
                value += conj(bra[local_index, j1, j2]) *
                         ket[local_index, source_j1, source_j2]
            end
        end
    end
    return value
end

"""
Lift one isolated magnetic Bloch band at every magnetic momentum to the
common six-orbital real-space grid.  These are full magnetic Bloch states,
not the zero-field band-projection weights used for visualization.
"""
function lift_magnetic_real_space_states(
    result::MagneticHWModel,
    band::Int;
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
    translation_start::Int=-div(result.l1, 2),
)
    1 <= band <= size(result.spectrum, 1) ||
        error("Magnetic band index lies outside the spectrum")
    model = result.model
    l1, l2, q = result.l1, result.l2, result.flux.q
    nlocal = local_dimension(model)
    fine_dimension1 = l1 * model.lg
    fine_dimension2 = l2 * q * model.lg
    translations = translation_start:(translation_start + l1 - 1)
    shift_maps = magnetic_shift_maps(
        result,
        translations;
        reciprocal_shift_periodic=reciprocal_shift_periodic,
    )

    states = Array{ComplexF32}(
        undef,
        nlocal,
        fine_dimension1,
        fine_dimension2,
        l1,
        l2,
    )
    plane_wave_amplitudes = zeros(ComplexF64, basis_dimension(model), 1)
    fine_momentum = zeros(ComplexF64, nlocal, fine_dimension1, fine_dimension2)
    real_space = similar(fine_momentum)

    for i2 in 1:l2, external_i1 in 1:l1
        ik = external_i1 + (i2 - 1) * l1
        coefficients = result.uort[:, :, ik] *
                       @view(result.eigvecs[:, band:band, ik])
        fill!(fine_momentum, 0)

        for ir_row in 1:q, i1_row in 1:l1
            TMBGCommonBasis.lifted_row_amplitudes!(
                plane_wave_amplitudes,
                result,
                coefficients,
                external_i1,
                i2,
                i1_row,
                ir_row,
                translations,
                shift_maps,
            )
            i2_full = (ir_row - 1) * l2 + i2
            for (g_index, g) in enumerate(model.Glist)
                fine_i1 = mod((i1_row - 1) + l1 * g[1], fine_dimension1) + 1
                fine_i2 = mod((i2_full - 1) + l2 * q * g[2], fine_dimension2) + 1
                rows = (nlocal * (g_index - 1) + 1):(nlocal * g_index)
                @views fine_momentum[:, fine_i1, fine_i2] .=
                    plane_wave_amplitudes[rows, 1]
            end
        end

        # Julia's two-dimensional ifft carries 1/(N1*N2).  Multiplication by
        # sqrt(N1*N2) makes the rectangular transform unitary.
        for local_index in 1:nlocal
            @views real_space[local_index, :, :] .=
                ifft(fine_momentum[local_index, :, :]) .*
                sqrt(fine_dimension1 * fine_dimension2)
        end
        normalization = sqrt(sum(abs2, real_space))
        normalization > 1e-12 || error("Lifted magnetic state has zero norm")
        @views states[:, :, :, external_i1, i2] .=
            ComplexF32.(real_space ./ normalization)
    end
    return states
end

function real_space_positions(
    result::MagneticHWModel,
    fine_dimension1::Int,
    fine_dimension2::Int,
)
    positions_x = zeros(Float64, fine_dimension1, fine_dimension2)
    positions_y = similar(positions_x)
    n1 = result.l1
    n2 = result.l2 * result.flux.q
    for j2 in 0:(fine_dimension2 - 1), j1 in 0:(fine_dimension1 - 1)
        centered_j1 = j1 < div(fine_dimension1 + 1, 2) ?
            j1 : j1 - fine_dimension1
        centered_j2 = j2 < div(fine_dimension2 + 1, 2) ?
            j2 : j2 - fine_dimension2
        position = (n1 / fine_dimension1) .* centered_j1 .* result.model.a1 .+
                   (n2 / fine_dimension2) .* centered_j2 .* result.model.a2
        positions_x[j1 + 1, j2 + 1] = position[1]
        positions_y[j1 + 1, j2 + 1] = position[2]
    end
    return positions_x, positions_y
end

function transport_phase(
    positions_x::Matrix{Float64},
    positions_y::Matrix{Float64},
    delta::Vector{Float64},
    reciprocal_wrap::Vector{Float64},
)
    effective_delta = delta .- reciprocal_wrap
    return @. exp(-1im * (
        effective_delta[1] * positions_x +
        effective_delta[2] * positions_y
    ))
end

function transported_overlap(
    bra::AbstractArray{<:Complex,3},
    ket::AbstractArray{<:Complex,3},
    phase::Matrix{ComplexF64},
)
    value = 0.0 + 0.0im
    for local_index in axes(bra, 1)
        @views value += sum(conj.(bra[local_index, :, :]) .*
                            phase .* ket[local_index, :, :])
    end
    return value
end

"""
Recover the Cartesian Fubini-Study metric from three gauge-invariant
projector distances measured along two, possibly nonorthogonal, momentum
steps.  `distance_difference` is the distance between the projectors at
`k + delta1` and `k + delta2`, so it probes the displacement
`delta2 - delta1`.

For a rank-`r` isolated subspace each distance is
`r - tr(P(k) P(k + delta))`.  The present magnetic-band calculation has
`r = 1`, for which this reduces to `1 - abs2(overlap)`.
"""
function cartesian_metric_from_projector_distances(
    delta1::AbstractVector{<:Real},
    delta2::AbstractVector{<:Real},
    distance1::Real,
    distance2::Real,
    distance_difference::Real,
)
    basis_step = hcat(delta1, delta2)
    abs(det(basis_step)) > 1e-15 || error("Degenerate momentum steps")
    cross = (distance1 + distance2 - distance_difference) / 2
    directional_metric = [distance1 cross; cross distance2]
    inverse_step = inv(basis_step)
    return transpose(inverse_step) * directional_metric * inverse_step
end

"""
Compute the gauge-invariant plaquette Berry curvature and the Fubini-Study
metric of an isolated magnetic hybrid-Wannier band.  The metric is first
obtained along the nonorthogonal reciprocal directions and then transformed
to Cartesian coordinates before taking Tr g.
"""
function compute_magnetic_geometry(
    result::MagneticHWModel,
    band::Int;
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
    translation_start::Int=-div(result.l1, 2),
    boundary_chern::Int=0,
)
    l1, l2, q = result.l1, result.l2, result.flux.q
    translations = translation_start:(translation_start + l1 - 1)
    shift_maps = magnetic_shift_maps(
        result,
        translations;
        reciprocal_shift_periodic=reciprocal_shift_periodic,
    )

    delta1 = result.model.g1 ./ l1
    delta2 = result.model.g2 ./ (l2 * q)
    basis_step = hcat(delta1, delta2)
    signed_area = det(basis_step)
    plaquette_area = abs(signed_area)
    abs(signed_area) > 1e-15 || error("Degenerate magnetic momentum mesh")

    link1 = zeros(ComplexF64, l1, l2)
    link2 = similar(link1)
    link_difference = similar(link1)

    current_row = lift_magnetic_momentum_row(
        result,
        band,
        1,
        translations,
        shift_maps,
    )
    for i2 in 1:l2
        next_row = lift_magnetic_momentum_row(
            result,
            band,
            i2 == l2 ? 1 : i2 + 1,
            translations,
            shift_maps,
        )
        for i1 in 1:l1
            j1 = mod1(i1 + 1, l1)
            current = @view current_row[:, :, :, i1]
            link1[i1, i2] = transported_overlap_momentum(
                current,
                @view(current_row[:, :, :, j1]),
                1,
                0,
            )
            link2[i1, i2] = transported_overlap_momentum(
                current,
                @view(next_row[:, :, :, i1]),
                0,
                1,
            )
            # Projector distance between k + delta1 and k + delta2.  This is
            # the finite-difference mixed term used by the LL implementation:
            # d(delta2-delta1) = d1 + d2 - 2 delta1' * g * delta2.
            link_difference[i1, i2] = transported_overlap_momentum(
                @view(current_row[:, :, :, j1]),
                @view(next_row[:, :, :, i1]),
                -1,
                1,
            )
        end
        current_row = next_row
        GC.gc(false)
    end

    # Optional explicit transition-function winding for alternative magnetic
    # gauges.  The production lifted-state convention already contains the
    # boundary sewing, so its default winding is zero.
    if boundary_chern != 0
        for i1 in 1:l1
            link2[i1, l2] *= exp(
                2im * pi * boundary_chern * (i1 - 1) / l1,
            )
        end
    end

    minimum(abs.(link1)) > 1e-10 || error("Vanishing x-link overlap")
    minimum(abs.(link2)) > 1e-10 || error("Vanishing y-link overlap")
    unit1 = link1 ./ abs.(link1)
    unit2 = link2 ./ abs.(link2)
    berry_flux = zeros(Float64, l1, l2)
    berry_curvature = similar(berry_flux)
    metric_trace = similar(berry_flux)
    trace_violation = similar(berry_flux)
    metric_xx = similar(berry_flux)
    metric_yy = similar(berry_flux)
    metric_xy = similar(berry_flux)

    for i2 in 1:l2, i1 in 1:l1
        j1 = mod1(i1 + 1, l1)
        j2 = mod1(i2 + 1, l2)
        loop = unit1[i1, i2] * unit2[j1, i2] *
               conj(unit1[i1, j2]) * conj(unit2[i1, i2])
        flux = angle(loop)
        berry_flux[i1, i2] = flux
        berry_curvature[i1, i2] = flux / signed_area

        distance1 = max(0.0, 1 - abs2(link1[i1, i2]))
        distance2 = max(0.0, 1 - abs2(link2[i1, i2]))
        distance_difference = max(
            0.0,
            1 - abs2(link_difference[i1, i2]),
        )
        cartesian_metric = cartesian_metric_from_projector_distances(
            delta1,
            delta2,
            distance1,
            distance2,
            distance_difference,
        )
        metric_xx[i1, i2] = cartesian_metric[1, 1]
        metric_yy[i1, i2] = cartesian_metric[2, 2]
        metric_xy[i1, i2] = cartesian_metric[1, 2]
        metric_trace[i1, i2] = tr(cartesian_metric)
        trace_violation[i1, i2] =
            metric_trace[i1, i2] - abs(berry_curvature[i1, i2])
    end

    chern = sum(berry_flux) / (2pi)
    eta = sum(trace_violation) * plaquette_area
    integrated_trace = sum(metric_trace) * plaquette_area
    integrated_abs_curvature = sum(abs, berry_curvature) * plaquette_area
    return MagneticGeometryResult(
        band,
        berry_flux,
        berry_curvature,
        metric_trace,
        trace_violation,
        metric_xx,
        metric_yy,
        metric_xy,
        abs.(link1),
        abs.(link2),
        chern,
        eta,
        integrated_trace,
        integrated_abs_curvature,
        plaquette_area,
    )
end

end
