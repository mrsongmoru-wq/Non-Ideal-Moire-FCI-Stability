module TMBGSymmetricGaugeProjection

using LinearAlgebra
using FFTW

using ..TMBGMagneticHW
using ..TMBGMagneticHW.TMBGHybridWannier.WangContinuum
using ..TMBGCommonBasis

export project_group_after_symmetric_gauge

function unitary_dft(dimension::Int)
    return [exp(2im * pi * row * column / dimension) / sqrt(dimension)
            for row in 0:(dimension - 1), column in 0:(dimension - 1)]
end

function symmetric_gauge_phase(
    result::MagneticHWModel,
    fine_dimension::Int,
    gauge_strength::Float64,
)
    l1 = result.l1
    lB = magnetic_length(result.model, result.flux.p, result.flux.q)
    phase = zeros(ComplexF64, fine_dimension, fine_dimension)
    for j2 in 0:(fine_dimension - 1), j1 in 0:(fine_dimension - 1)
        centered_j1 = j1 < div(fine_dimension + 1, 2) ? j1 : j1 - fine_dimension
        centered_j2 = j2 < div(fine_dimension + 1, 2) ? j2 : j2 - fine_dimension
        position = (l1 / fine_dimension) .* (
            centered_j1 .* result.model.a1 .+
            centered_j2 .* result.model.a2
        )
        exponent = -result.flux.signB * gauge_strength *
                   position[1] * position[2] / (2lB^2)
        phase[j1 + 1, j2 + 1] = exp(1im * exponent)
    end
    return phase
end

function project_group_after_symmetric_gauge(
    result::MagneticHWModel,
    selected_bands::Vector{Int};
    target_local_band::Int=div(size(result.hw.WLS, 2), 2) + 1,
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
    gauge_strength::Float64=1.0,
    transform_method::Symbol=:fft,
    translation_start::Int=-div(result.l1, 2),
    metric_tolerance::Float64=1e-12,
)
    model = result.model
    nactive = size(result.hw.WLS, 2)
    nlocal = local_dimension(model)
    l1, l2, q = result.l1, result.l2, result.flux.q
    nselected = length(selected_bands)
    nk = l1 * l2
    translations = translation_start:(translation_start + l1 - 1)
    fine_dimension = l1 * model.lg

    isempty(selected_bands) && error("At least one magnetic band must be selected")
    all((1 .<= selected_bands) .& (selected_bands .<= size(result.spectrum, 1))) ||
        error("Selected magnetic band lies outside the spectrum")
    1 <= target_local_band <= nactive ||
        error("Target zero-field band lies outside the active space")
    transform_method in (:fft, :direct) ||
        error("transform_method must be :fft or :direct")

    p_signed, nq = result.flux.signB * result.flux.p, result.nq
    shifts = Set{Tuple{Int,Int}}()
    for ir_row in 1:q, i1_row in 1:l1, t in translations
        i1_unwrapped = i1_row - t * p_signed * nq
        ir_unwrapped = ir_row - t * p_signed
        push!(shifts, (
            div(i1_unwrapped - mod1(i1_unwrapped, l1), l1),
            div(ir_unwrapped - mod1(ir_unwrapped, q), q),
        ))
    end
    shift_maps = Dict(shift => TMBGCommonBasis.reciprocal_shift_map(
        model,
        shift...;
        periodic=reciprocal_shift_periodic,
    ) for shift in shifts)

    fourier = transform_method == :direct ? unitary_dft(fine_dimension) : nothing
    gauge_phase = symmetric_gauge_phase(result, fine_dimension, gauge_strength)
    band_weights = zeros(Float64, nactive, l1, l1)
    common_weights = zeros(Float64, l1, l1)
    raw_metric_errors = zeros(Float64, nk)
    raw_metric_eigenvalue_bounds = zeros(Float64, nk, 2)
    plane_wave_amplitudes = zeros(
        ComplexF64,
        basis_dimension(model),
        nselected,
    )
    fine_momentum = zeros(
        ComplexF64,
        nlocal,
        fine_dimension,
        fine_dimension,
        nselected,
    )
    transformed_momentum = similar(fine_momentum)
    bloch_amplitudes = zeros(ComplexF64, basis_dimension(model), nselected)
    energy_amplitudes = zeros(ComplexF64, nactive, nselected)

    for i2 in 1:l2, external_i1 in 1:l1
        ik = external_i1 + (i2 - 1) * l1
        coefficients = result.uort[:, :, ik] *
                       result.eigvecs[:, selected_bands, ik]
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
                fine_i1 = mod((i1_row - 1) + l1 * g[1], fine_dimension) + 1
                fine_i2 = mod((i2_full - 1) + l1 * g[2], fine_dimension) + 1
                rows = (nlocal * (g_index - 1) + 1):(nlocal * g_index)
                fine_momentum[:, fine_i1, fine_i2, :] .=
                    @view(plane_wave_amplitudes[rows, :])
            end
        end

        flattened = reshape(fine_momentum, :, nselected)
        metric = Hermitian(flattened' * flattened)
        metric_values = eigvals(metric)
        minimum(metric_values) > metric_tolerance ||
            error("Common plane-wave image is rank deficient at magnetic k index $ik")
        raw_metric_errors[ik] = norm(Matrix(metric) - I) / sqrt(nselected)
        raw_metric_eigenvalue_bounds[ik, :] .= extrema(metric_values)
        correction = TMBGCommonBasis.inverse_sqrt_metric(
            metric;
            tolerance=metric_tolerance,
        )
        flattened .= flattened * correction

        for selected_index in 1:nselected, local_index in 1:nlocal
            momentum_slice = @view fine_momentum[local_index, :, :, selected_index]
            if transform_method == :fft
                real_space = ifft(momentum_slice) .* fine_dimension
                transformed_momentum[local_index, :, :, selected_index] .=
                    fft(real_space .* gauge_phase) ./ fine_dimension
            else
                real_space = fourier * momentum_slice * transpose(fourier)
                transformed_momentum[local_index, :, :, selected_index] .=
                    fourier' * (real_space .* gauge_phase) * conj(fourier)
            end
        end

        for i2_full in 1:l1, i1_row in 1:l1
            fill!(bloch_amplitudes, 0)
            for (g_index, g) in enumerate(model.Glist)
                fine_i1 = mod((i1_row - 1) + l1 * g[1], fine_dimension) + 1
                fine_i2 = mod((i2_full - 1) + l1 * g[2], fine_dimension) + 1
                rows = (nlocal * (g_index - 1) + 1):(nlocal * g_index)
                bloch_amplitudes[rows, :] .= reshape(
                    @view(transformed_momentum[:, fine_i1, fine_i2, :]),
                    nlocal,
                    nselected,
                )
            end
            _, energy = TMBGMagneticHW.TMBGHybridWannier.energy_eigenbasis(
                result.hw,
                i1_row,
                i2_full,
            )
            mul!(energy_amplitudes, energy', bloch_amplitudes)
            band_weights[:, i1_row, i2_full] .+=
                vec(sum(abs2.(energy_amplitudes); dims=2))
            common_weights[i1_row, i2_full] += sum(abs2, bloch_amplitudes)
        end
    end

    band_weights ./= nk * nselected
    common_weights ./= nk * nselected
    target_weights = Matrix(@view band_weights[target_local_band, :, :])
    active_weights = dropdims(sum(band_weights; dims=1); dims=1)
    source_weights = vec(dropdims(sum(band_weights; dims=(2, 3)); dims=(2, 3)))
    target_fraction = sum(target_weights) /
        max(sum(active_weights), eps(Float64))
    return TMBGCommonBasis.CommonBasisProjection(
        selected_bands,
        band_weights,
        target_weights,
        active_weights,
        common_weights,
        source_weights,
        target_fraction,
        raw_metric_errors,
        raw_metric_eigenvalue_bounds,
    )
end

end
