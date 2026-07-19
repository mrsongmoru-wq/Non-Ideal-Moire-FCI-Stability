module TMBGIdealComponent

using LinearAlgebra
using Statistics

using ..TMBGMagneticHW
using ..TMBGMagneticHW.TMBGHybridWannier.WangContinuum
using ..TMBGProjection
using ..TMBGCommonBasis
using ..TMBGMagneticGeometry
using ..TMBGIntrinsicIdealGeometry

export IdealComponentResult,
       select_isolated_target_band,
       compute_ideal_component,
       write_trace_condition_results,
       write_overlap_diagnostics

"""
Geometry of one target magnetic subband after a multiband correction frame
has removed the nonunitary distortion of the truncated hybrid-Wannier basis.

The target band is always the first column of `correction_bands`.  The
remaining columns are nearby magnetic subbands used only to define a smooth
covariant frame; they are not combined into the target projector.
"""
struct IdealComponentResult
    target_band::Int
    correction_bands::Vector{Int}
    berry_flux::Matrix{Float64}
    berry_curvature::Matrix{Float64}
    ideal_metric_xx::Matrix{Float64}
    ideal_metric_yy::Matrix{Float64}
    ideal_metric_xy::Matrix{Float64}
    residual_metric_xx::Matrix{Float64}
    residual_metric_yy::Matrix{Float64}
    residual_metric_xy::Matrix{Float64}
    residual_metric_trace::Matrix{Float64}
    intrinsic_metric_xx::Matrix{Float64}
    intrinsic_metric_yy::Matrix{Float64}
    intrinsic_metric_xy::Matrix{Float64}
    intrinsic_metric_trace::Matrix{Float64}
    trace_excess::Matrix{Float64}
    raw_link1_abs::Matrix{Float64}
    raw_link2_abs::Matrix{Float64}
    covariant_link1_abs::Matrix{Float64}
    covariant_link2_abs::Matrix{Float64}
    frame_singular1::Matrix{Float64}
    frame_singular2::Matrix{Float64}
    frame_singular_difference::Matrix{Float64}
    chern::Float64
    frame_chern::Float64
    ideal_integrated_trace::Float64
    residual_integrated_trace::Float64
    total_integrated_trace::Float64
    eta::Float64
    ideal_fraction::Float64
    minimum_frame_link_singular_value::Float64
    plaquette_area::Float64
end

function select_isolated_target_band(result::MagneticHWModel)
    _, upper = select_target_upper_group(
        result.spectrum,
        result.flux.q;
        target_chern=-2,
        upper_reference_energy=25.047077583,
    )
    length(upper) == 1 ||
        error("Expected one isolated upper target branch, got $upper")
    return only(upper)
end

function target_preserving_orthonormalize!(states::Matrix{ComplexF64})
    for column in axes(states, 2)
        for previous in 1:(column - 1)
            @views states[:, column] .-=
                dot(states[:, previous], states[:, column]) .* states[:, previous]
        end
        normalization = norm(@view states[:, column])
        normalization > 1e-10 || error("Correction frame became rank deficient")
        @views states[:, column] ./= normalization

        # Reorthogonalization controls loss of precision for nearby remote
        # bands while the first (target) column remains unchanged.
        for previous in 1:(column - 1)
            @views states[:, column] .-=
                dot(states[:, previous], states[:, column]) .* states[:, previous]
        end
        normalization = norm(@view states[:, column])
        normalization > 1e-10 || error("Correction frame became rank deficient")
        @views states[:, column] ./= normalization
    end
    return states
end

function lift_frame(
    result::MagneticHWModel,
    ordered_bands::Vector{Int},
    external_i1::Int,
    i2::Int,
    translations::UnitRange{Int},
    shift_maps::Dict{Tuple{Int,Int},Vector{Int}},
)
    model = result.model
    nlocal = local_dimension(model)
    dimension1 = result.l1 * model.lg
    dimension2 = result.l2 * result.flux.q * model.lg
    ik = external_i1 + (i2 - 1) * result.l1
    coefficients = result.uort[:, :, ik] *
                   @view(result.eigvecs[:, ordered_bands, ik])
    nstates = length(ordered_bands)
    plane_wave = zeros(ComplexF64, basis_dimension(model), nstates)
    momentum = zeros(ComplexF64, nlocal, dimension1, dimension2, nstates)

    for ir_row in 1:result.flux.q, i1_row in 1:result.l1
        TMBGCommonBasis.lifted_row_amplitudes!(
            plane_wave,
            result,
            coefficients,
            external_i1,
            i2,
            i1_row,
            ir_row,
            translations,
            shift_maps,
        )
        i2_full = (ir_row - 1) * result.l2 + i2
        for (g_index, g) in enumerate(model.Glist)
            fine_i1 = mod((i1_row - 1) + result.l1 * g[1], dimension1) + 1
            fine_i2 = mod(
                (i2_full - 1) + result.l2 * result.flux.q * g[2],
                dimension2,
            ) + 1
            rows = (nlocal * (g_index - 1) + 1):(nlocal * g_index)
            @views momentum[:, fine_i1, fine_i2, :] .= plane_wave[rows, :]
        end
    end
    states = reshape(momentum, :, nstates)
    return target_preserving_orthonormalize!(states)
end

function momentum_permutation(
    nlocal::Int,
    dimension1::Int,
    dimension2::Int,
    shift1::Int,
    shift2::Int,
)
    permutation = Vector{Int}(undef, nlocal * dimension1 * dimension2)
    destination = 1
    for j2 in 1:dimension2, j1 in 1:dimension1, local_index in 1:nlocal
        source_j1 = mod1(j1 + shift1, dimension1)
        source_j2 = mod1(j2 + shift2, dimension2)
        permutation[destination] =
            local_index + nlocal * (source_j1 - 1) +
            nlocal * dimension1 * (source_j2 - 1)
        destination += 1
    end
    return permutation
end

function frame_link(
    bra::Matrix{ComplexF64},
    ket::Matrix{ComplexF64},
    permutation::Vector{Int},
)
    overlap = bra' * ket[permutation, :]
    transport, singular_values = polar_frame_transport(overlap)
    return overlap[1, 1], transport[1, 1], minimum(singular_values), det(transport)
end

function link_geometry(
    link1::Matrix{ComplexF64},
    link2::Matrix{ComplexF64},
    link_difference::Matrix{ComplexF64},
    delta1::Vector{Float64},
    delta2::Vector{Float64},
)
    l1, l2 = size(link1)
    signed_area = det(hcat(delta1, delta2))
    abs(signed_area) > 1e-15 || error("Degenerate magnetic momentum mesh")
    plaquette_area = abs(signed_area)
    unit1 = link1 ./ abs.(link1)
    unit2 = link2 ./ abs.(link2)
    flux = zeros(Float64, l1, l2)
    metric_xx = similar(flux)
    metric_yy = similar(flux)
    metric_xy = similar(flux)
    metric_trace = similar(flux)

    for i2 in 1:l2, i1 in 1:l1
        loop = unit1[i1, i2] * unit2[mod1(i1 + 1, l1), i2] *
               conj(unit1[i1, mod1(i2 + 1, l2)]) * conj(unit2[i1, i2])
        flux[i1, i2] = angle(loop)
        metric = cartesian_metric_from_projector_distances(
            delta1,
            delta2,
            max(0.0, 1 - abs2(link1[i1, i2])),
            max(0.0, 1 - abs2(link2[i1, i2])),
            max(0.0, 1 - abs2(link_difference[i1, i2])),
        )
        metric_xx[i1, i2] = metric[1, 1]
        metric_yy[i1, i2] = metric[2, 2]
        metric_xy[i1, i2] = metric[1, 2]
        metric_trace[i1, i2] = tr(metric)
    end
    curvature = flux ./ signed_area
    return (
        flux=flux,
        curvature=curvature,
        metric_xx=metric_xx,
        metric_yy=metric_yy,
        metric_xy=metric_xy,
        metric_trace=metric_trace,
        chern=sum(flux) / (2pi),
        plaquette_area=plaquette_area,
    )
end

"""
Compute the intrinsic ideal component of a magnetic hybrid-Wannier subband.

At each momentum, the target state is completed by nearby magnetic subbands.
The polar part of the inter-frame overlap supplies a gauge-covariant target
connection free of the nonunitary distortion of the finite hWF frame.  Berry
curvature follows from its Wilson plaquette.  Three projector distances along
`delta1`, `delta2`, and `delta2-delta1` determine the Cartesian residual
metric, including the angle between reciprocal-lattice directions.

For the C3-symmetric tMBG model, the local Kähler component is uniquely
`g_ideal = |Omega| I/2`.  The reported intrinsic metric is the sum of this
ideal component and the covariant residual.  Consequently
`eta = integral Tr(g_residual)` and the ideal fraction is
`integral |Omega| / integral Tr(g_intrinsic)`.
"""
function compute_ideal_component(
    result::MagneticHWModel;
    target_band::Int=select_isolated_target_band(result),
    correction_bands::Vector{Int}=[target_band - 1],
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
    translation_start::Int=-div(result.l1, 2),
)
    requested = unique(vcat(target_band, correction_bands))
    ordered_bands = vcat(target_band, sort(setdiff(requested, [target_band])))
    length(ordered_bands) >= 2 ||
        error("Use at least one correction band for the multiband frame")
    all(1 .<= ordered_bands .<= size(result.spectrum, 1)) ||
        error("A correction band lies outside the magnetic spectrum")

    translations = translation_start:(translation_start + result.l1 - 1)
    shift_maps = TMBGMagneticGeometry.magnetic_shift_maps(
        result,
        translations;
        reciprocal_shift_periodic=reciprocal_shift_periodic,
    )
    nlocal = local_dimension(result.model)
    dimension1 = result.l1 * result.model.lg
    dimension2 = result.l2 * result.flux.q * result.model.lg
    permutation1 = momentum_permutation(nlocal, dimension1, dimension2, 1, 0)
    permutation2 = momentum_permutation(nlocal, dimension1, dimension2, 0, 1)
    permutation_difference =
        momentum_permutation(nlocal, dimension1, dimension2, -1, 1)

    l1, l2, q = result.l1, result.l2, result.flux.q
    raw1 = zeros(ComplexF64, l1, l2)
    raw2 = similar(raw1)
    raw_difference = similar(raw1)
    covariant1 = similar(raw1)
    covariant2 = similar(raw1)
    covariant_difference = similar(raw1)
    frame_determinant1 = similar(raw1)
    frame_determinant2 = similar(raw1)
    singular1_map = zeros(Float64, l1, l2)
    singular2_map = similar(singular1_map)
    singular_difference_map = similar(singular1_map)

    for i2 in 1:l2
        next_i2 = mod1(i2 + 1, l2)
        current = lift_frame(
            result,
            ordered_bands,
            1,
            i2,
            translations,
            shift_maps,
        )
        next_row = lift_frame(
            result,
            ordered_bands,
            1,
            next_i2,
            translations,
            shift_maps,
        )
        for i1 in 1:l1
            current_x = lift_frame(
                result,
                ordered_bands,
                mod1(i1 + 1, l1),
                i2,
                translations,
                shift_maps,
            )
            raw1[i1, i2], covariant1[i1, i2], singular1_map[i1, i2],
                frame_determinant1[i1, i2] =
                frame_link(current, current_x, permutation1)
            raw2[i1, i2], covariant2[i1, i2], singular2_map[i1, i2],
                frame_determinant2[i1, i2] =
                frame_link(current, next_row, permutation2)
            raw_difference[i1, i2], covariant_difference[i1, i2],
                singular_difference_map[i1, i2], _ =
                frame_link(current_x, next_row, permutation_difference)
            current = current_x
            if i1 < l1
                next_row = lift_frame(
                    result,
                    ordered_bands,
                    i1 + 1,
                    next_i2,
                    translations,
                    shift_maps,
                )
            end
            GC.gc(false)
        end
    end

    delta1 = result.model.g1 ./ l1
    delta2 = result.model.g2 ./ (l2 * q)
    covariant = link_geometry(
        covariant1,
        covariant2,
        covariant_difference,
        delta1,
        delta2,
    )

    ideal_xx = zeros(Float64, l1, l2)
    ideal_yy = similar(ideal_xx)
    ideal_xy = similar(ideal_xx)
    for index in eachindex(covariant.curvature)
        ideal_metric = c3_kahler_metric(covariant.curvature[index])
        ideal_xx[index] = ideal_metric[1, 1]
        ideal_yy[index] = ideal_metric[2, 2]
        ideal_xy[index] = ideal_metric[1, 2]
    end
    intrinsic_xx = ideal_xx .+ covariant.metric_xx
    intrinsic_yy = ideal_yy .+ covariant.metric_yy
    intrinsic_xy = ideal_xy .+ covariant.metric_xy
    intrinsic_trace = abs.(covariant.curvature) .+ covariant.metric_trace
    trace_excess = copy(covariant.metric_trace)

    ideal_trace = sum(abs, covariant.curvature) * covariant.plaquette_area
    residual_trace = sum(trace_excess) * covariant.plaquette_area
    total_trace = sum(intrinsic_trace) * covariant.plaquette_area

    frame_unit1 = frame_determinant1 ./ abs.(frame_determinant1)
    frame_unit2 = frame_determinant2 ./ abs.(frame_determinant2)
    frame_flux = 0.0
    for i2 in 1:l2, i1 in 1:l1
        frame_flux += angle(
            frame_unit1[i1, i2] * frame_unit2[mod1(i1 + 1, l1), i2] *
            conj(frame_unit1[i1, mod1(i2 + 1, l2)]) *
            conj(frame_unit2[i1, i2]),
        )
    end
    minimum_singular = min(
        minimum(singular1_map),
        minimum(singular2_map),
        minimum(singular_difference_map),
    )

    return IdealComponentResult(
        target_band,
        ordered_bands,
        covariant.flux,
        covariant.curvature,
        ideal_xx,
        ideal_yy,
        ideal_xy,
        covariant.metric_xx,
        covariant.metric_yy,
        covariant.metric_xy,
        covariant.metric_trace,
        intrinsic_xx,
        intrinsic_yy,
        intrinsic_xy,
        intrinsic_trace,
        trace_excess,
        abs.(raw1),
        abs.(raw2),
        abs.(covariant1),
        abs.(covariant2),
        singular1_map,
        singular2_map,
        singular_difference_map,
        covariant.chern,
        frame_flux / (2pi),
        ideal_trace,
        residual_trace,
        total_trace,
        residual_trace,
        ideal_trace / total_trace,
        minimum_singular,
        covariant.plaquette_area,
    )
end

function write_trace_condition_results(
    geometry::IdealComponentResult,
    result::MagneticHWModel,
    output_directory::String,
)
    mkpath(output_directory)
    summary_path = joinpath(output_directory, "trace_condition_summary.txt")
    csv_path = joinpath(output_directory, "quantum_geometry.csv")
    open(summary_path, "w") do io
        println(io, "method=multiband_polar_frame_plus_C3_Kahler_component")
        println(io, "p=$(result.flux.p)")
        println(io, "q=$(result.flux.q)")
        println(io, "target_band=$(geometry.target_band)")
        println(io, "correction_frame=$(join(geometry.correction_bands, ','))")
        println(io, "frame_rank=$(length(geometry.correction_bands))")
        println(io, "chern=$(geometry.chern)")
        println(io, "frame_chern=$(geometry.frame_chern)")
        println(io, "minimum_frame_link_singular_value=$(geometry.minimum_frame_link_singular_value)")
        println(io, "ideal_integrated_trace=$(geometry.ideal_integrated_trace)")
        println(io, "residual_integrated_trace=$(geometry.residual_integrated_trace)")
        println(io, "total_integrated_trace=$(geometry.total_integrated_trace)")
        println(io, "eta=$(geometry.eta)")
        println(io, "ideal_fraction=$(geometry.ideal_fraction)")
        println(io, "berry_curvature_min=$(minimum(geometry.berry_curvature))")
        println(io, "berry_curvature_max=$(maximum(geometry.berry_curvature))")
    end

    l1, l2 = size(geometry.berry_curvature)
    open(csv_path, "w") do io
        println(
            io,
            "k1,k2,berry_curvature,ideal_g_xx,ideal_g_yy,ideal_g_xy," *
            "residual_g_xx,residual_g_yy,residual_g_xy,residual_trace," *
            "intrinsic_g_xx,intrinsic_g_yy,intrinsic_g_xy,intrinsic_trace," *
            "trace_excess",
        )
        for i2 in 1:l2, i1 in 1:l1
            println(io, join((
                (i1 - 1) / l1,
                (i2 - 1) / (l2 * result.flux.q),
                geometry.berry_curvature[i1, i2],
                geometry.ideal_metric_xx[i1, i2],
                geometry.ideal_metric_yy[i1, i2],
                geometry.ideal_metric_xy[i1, i2],
                geometry.residual_metric_xx[i1, i2],
                geometry.residual_metric_yy[i1, i2],
                geometry.residual_metric_xy[i1, i2],
                geometry.residual_metric_trace[i1, i2],
                geometry.intrinsic_metric_xx[i1, i2],
                geometry.intrinsic_metric_yy[i1, i2],
                geometry.intrinsic_metric_xy[i1, i2],
                geometry.intrinsic_metric_trace[i1, i2],
                geometry.trace_excess[i1, i2],
            ), ','))
        end
    end
    return summary_path, csv_path
end

function write_distribution_line(io, label::String, values::AbstractArray{<:Real})
    flattened = vec(values)
    println(io, "$(label)_min=$(minimum(flattened))")
    println(io, "$(label)_q05=$(quantile(flattened, 0.05))")
    println(io, "$(label)_median=$(median(flattened))")
    println(io, "$(label)_mean=$(mean(flattened))")
    println(io, "$(label)_q95=$(quantile(flattened, 0.95))")
    println(io, "$(label)_max=$(maximum(flattened))")
    println(io, "$(label)_std=$(std(flattened; corrected=false))")
end

"""
Write every target-link modulus and correction-frame singular value on the
magnetic Brillouin-zone mesh.  These data expose the full overlap distribution
behind the ideal-component calculation rather than only its extrema.
"""
function write_overlap_diagnostics(
    geometry::IdealComponentResult,
    result::MagneticHWModel,
    output_directory::String,
)
    mkpath(output_directory)
    csv_path = joinpath(output_directory, "overlap_diagnostics.csv")
    summary_path = joinpath(output_directory, "overlap_diagnostics_summary.txt")
    l1, l2 = size(geometry.raw_link1_abs)
    open(csv_path, "w") do io
        println(
            io,
            "k1,k2,raw_link1_abs,raw_link2_abs,covariant_link1_abs," *
            "covariant_link2_abs,frame_singular1,frame_singular2," *
            "frame_singular_difference",
        )
        for i2 in 1:l2, i1 in 1:l1
            println(io, join((
                (i1 - 1) / l1,
                (i2 - 1) / (l2 * result.flux.q),
                geometry.raw_link1_abs[i1, i2],
                geometry.raw_link2_abs[i1, i2],
                geometry.covariant_link1_abs[i1, i2],
                geometry.covariant_link2_abs[i1, i2],
                geometry.frame_singular1[i1, i2],
                geometry.frame_singular2[i1, i2],
                geometry.frame_singular_difference[i1, i2],
            ), ','))
        end
    end
    open(summary_path, "w") do io
        println(io, "target_band=$(geometry.target_band)")
        println(io, "correction_frame=$(join(geometry.correction_bands, ','))")
        for (label, values) in (
            ("raw_link1_abs", geometry.raw_link1_abs),
            ("raw_link2_abs", geometry.raw_link2_abs),
            ("covariant_link1_abs", geometry.covariant_link1_abs),
            ("covariant_link2_abs", geometry.covariant_link2_abs),
            ("frame_singular1", geometry.frame_singular1),
            ("frame_singular2", geometry.frame_singular2),
            ("frame_singular_difference", geometry.frame_singular_difference),
        )
            write_distribution_line(io, label, values)
        end
    end
    return summary_path, csv_path
end

end
