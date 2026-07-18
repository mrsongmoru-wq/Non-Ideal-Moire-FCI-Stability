module TMBGFigure4Projection

using DelimitedFiles
using LinearAlgebra
using Serialization
using Statistics

using ..TMBGMagneticHW
using ..TMBGMagneticHW.TMBGHybridWannier.WangContinuum
using ..TMBGProjection
using ..TMBGCommonBasis
using ..TMBGIdealComponent
using ..TMBGSymmetricGaugeProjection

export write_figure4_projection

function distribute_segments(nodes, model, total_points::Int)
    distances = Float64[]
    for i in 1:(length(nodes) - 1)
        delta = nodes[i + 1] .- nodes[i]
        cartesian = delta[1] .* model.g1 .+ delta[2] .* model.g2
        push!(distances, norm(cartesian))
    end
    intervals = max.(
        1,
        round.(Int, distances ./ sum(distances) .* (total_points - 1)),
    )
    while sum(intervals) != total_points - 1
        if sum(intervals) < total_points - 1
            intervals[argmax(distances ./ intervals)] += 1
        else
            candidates = findall(>(1), intervals)
            isempty(candidates) && error("Cannot reduce path intervals")
            index = candidates[argmin(distances[candidates] ./ intervals[candidates])]
            intervals[index] -= 1
        end
    end
    return intervals
end

function build_path(nodes, model; total_points::Int=241)
    intervals = distribute_segments(nodes, model, total_points)
    points = Vector{Vector{Float64}}()
    cumulative_distance = Float64[]
    node_indices = Int[1]
    distance = 0.0
    for segment in eachindex(intervals)
        count = intervals[segment]
        for step in 0:(count - 1)
            fraction = step / count
            point =
                (1 - fraction) .* nodes[segment] .+ fraction .* nodes[segment + 1]
            if !isempty(points)
                delta = point .- points[end]
                distance += norm(delta[1] .* model.g1 .+ delta[2] .* model.g2)
            end
            push!(points, point)
            push!(cumulative_distance, distance)
        end
        push!(node_indices, length(points) + 1)
    end
    final_point = copy(nodes[end])
    delta = final_point .- points[end]
    distance += norm(delta[1] .* model.g1 .+ delta[2] .* model.g2)
    push!(points, final_point)
    push!(cumulative_distance, distance)
    node_indices[end] = length(points)
    return points, cumulative_distance, node_indices
end

function nearest_grid_index(fraction::Real, dimension::Int)
    return mod(round(Int, fraction * dimension), dimension) + 1
end

function c3_orbit_average(weights::Matrix{Float64}, index::Tuple{Int,Int})
    dimension = size(weights, 1)
    i1, i2 = index[1] - 1, index[2] - 1
    total = 0.0
    for _ in 1:3
        total += weights[i1 + 1, i2 + 1]
        i1, i2 = mod(-i2, dimension), mod(i1 - i2, dimension)
    end
    return total / 3
end

function write_projection_grid(
    projection::CommonBasisProjection,
    result::MagneticHWModel,
    target_local::Int,
    outdir::String,
)
    l1 = result.l1
    target_sum = sum(projection.target_weights)
    target_max = maximum(projection.target_weights)
    rows = zeros(Float64, l1 * l1, 10)
    row = 1
    for i2 in 1:l1, i1 in 1:l1
        k1 = (i1 - 1) / l1
        k2 = (i2 - 1) / l1
        kcart = k1 .* result.model.g1 .+ k2 .* result.model.g2
        weight = projection.target_weights[i1, i2]
        rows[row, :] .= [
            k1,
            k2,
            kcart[1],
            kcart[2],
            result.hw.eigvals[target_local, i1, i2],
            weight,
            weight / max(target_sum, eps(Float64)),
            weight / max(target_max, eps(Float64)),
            projection.active_weights[i1, i2],
            projection.common_weights[i1, i2],
        ]
        row += 1
    end
    grid_path = joinpath(outdir, "upper_component_grid.csv")
    open(grid_path, "w") do io
        println(
            io,
            "k1_wang,k2_wang,kx,ky,zero_field_energy_meV,target_weight," *
            "target_probability,target_relative,active_weight,common_weight",
        )
        writedlm(io, rows, ',')
    end

    source_rows = hcat(
        collect(1:length(projection.source_weights)),
        result.hw.band_indices,
        projection.source_weights,
        projection.source_weights ./
        max(sum(projection.source_weights), eps(Float64)),
    )
    source_path = joinpath(outdir, "upper_zero_field_band_sources.csv")
    open(source_path, "w") do io
        println(
            io,
            "local_band,global_band,weight," *
            "fraction_within_retained_zero_field_bands",
        )
        writedlm(io, source_rows, ',')
    end
    return grid_path, source_path
end

function write_zero_field_path(
    result::MagneticHWModel,
    target_local::Int,
    outdir::String;
    total_points::Int=241,
)
    model = result.model
    labels = ["K", "Gamma", "M", "Kprime"]
    nodes = [
        [1 / 3, 2 / 3],
        [2 / 3, 1 / 3],
        [1 / 6, 1 / 3],
        [0.0, 0.0],
    ]
    points, distances, node_indices =
        build_path(nodes, model; total_points=total_points)
    active_bands = result.hw.band_indices
    energies = zeros(Float64, length(points), length(active_bands))
    for (index, reduced) in enumerate(points)
        momentum = reduced[1] .* model.g1 .+ reduced[2] .* model.g2
        spectrum = eigvals(Hermitian(WangContinuum.hamiltonian(model, momentum)))
        energies[index, :] .= spectrum[active_bands]
    end

    kpoint_path = joinpath(outdir, "zero_field_path_kpoints.csv")
    path_rows = hcat(
        collect(0:(length(points) - 1)),
        reduce(hcat, points)',
        distances,
    )
    open(kpoint_path, "w") do io
        println(io, "path_index,k1_wang,k2_wang,k_distance")
        writedlm(io, path_rows, ',')
    end
    energy_path = joinpath(outdir, "zero_field_path_energies.csv")
    open(energy_path, "w") do io
        println(io, join(["band_$(band)_meV" for band in active_bands], ','))
        writedlm(io, energies, ',')
    end
    metadata_path = joinpath(outdir, "zero_field_path_metadata.txt")
    open(metadata_path, "w") do io
        println(io, "basis=Wang_tMBG_six_orbital_plane_wave")
        println(io, "labels=$(join(labels, ','))")
        println(io, "node_indices_zero_based=$(join(node_indices .- 1, ','))")
        println(
            io,
            "node_coordinates_wang=$(join([join(node, ':') for node in nodes], ','))",
        )
        println(io, "active_global_bands=$(join(active_bands, ','))")
        println(io, "target_local_band=$target_local")
        println(io, "target_global_band=$(active_bands[target_local])")
        println(io, "g1=$(join(model.g1, ','))")
        println(io, "g2=$(join(model.g2, ','))")
        println(io, "num_path_points=$(length(points))")
    end
    return kpoint_path, energy_path, metadata_path
end

"""
Reproduce the momentum-resolved ideal-component projection used in Fig. 4.

The isolated weak-field magnetic branch is lifted to the common six-orbital
plane-wave Hilbert space, transformed from the hWF Landau-gauge convention to
symmetric gauge, and projected onto the true zero-field energy eigenbasis.
The exported map is the target-band weight normalized by its maximum; the
91.9% quantity is the integrated target weight within the central zero-field
pair.
"""
function write_figure4_projection(
    result::MagneticHWModel,
    outdir::String;
    target_band::Int=select_paper_target_band(result),
    total_path_points::Int=241,
)
    mkpath(outdir)
    target_local = div(size(result.hw.WLS, 2), 2) + 1
    projection = project_group_after_symmetric_gauge(
        result,
        [target_band];
        target_local_band=target_local,
        reciprocal_shift_periodic=true,
        gauge_strength=1.0,
        transform_method=:fft,
        translation_start=-div(result.l1, 2),
    )
    serialize(joinpath(outdir, "figure4_projection.jls"), projection)
    grid_path, source_path =
        write_projection_grid(projection, result, target_local, outdir)
    kpoint_path, energy_path, path_metadata = write_zero_field_path(
        result,
        target_local,
        outdir;
        total_points=total_path_points,
    )

    target_max = maximum(projection.target_weights)
    gamma_index = (
        nearest_grid_index(2 / 3, result.l1),
        nearest_grid_index(1 / 3, result.l1),
    )
    k_index = (
        nearest_grid_index(1 / 3, result.l1),
        nearest_grid_index(2 / 3, result.l1),
    )
    m_index = (
        nearest_grid_index(1 / 6, result.l1),
        nearest_grid_index(1 / 3, result.l1),
    )
    kprime_index = (1, 1)
    central_pair = (target_local - 1):target_local
    central_pair_fraction = projection.source_weights[target_local] /
                            sum(projection.source_weights[central_pair])
    metadata_path = joinpath(outdir, "metadata.txt")
    open(metadata_path, "w") do io
        println(io, "method=common_six_orbital_plane_wave_projection")
        println(io, "gauge_transform=exp(-i signB x y / (2 lB^2))")
        println(io, "selected_magnetic_bands=$target_band")
        println(io, "target_local_band=$target_local")
        println(io, "target_global_band=$(result.hw.band_indices[target_local])")
        println(io, "common_norm=$(sum(projection.common_weights))")
        println(io, "retained_zero_field_active_weight=$(sum(projection.active_weights))")
        println(
            io,
            "target_fraction_within_retained_active=" *
            "$(projection.target_fraction_within_active)",
        )
        println(io, "target_fraction_within_central_pair=$central_pair_fraction")
        println(io, "target_c3_relative_error=$(c3_relative_error(projection.target_weights))")
        println(io, "raw_metric_error_max=$(maximum(projection.raw_metric_errors))")
        for (label, index) in (
            ("Gamma", gamma_index),
            ("K", k_index),
            ("M", m_index),
            ("Kprime", kprime_index),
        )
            weight = projection.target_weights[index...]
            orbit_average = c3_orbit_average(projection.target_weights, index)
            println(io, "$(label)_target_weight=$weight")
            println(io, "$(label)_c3_orbit_average=$orbit_average")
            println(
                io,
                "$(label)_c3_orbit_average_relative_to_max=" *
                "$(orbit_average / target_max)",
            )
        end
        println(io, "maximum_target_weight=$target_max")
    end
    return (
        projection=projection,
        grid=grid_path,
        sources=source_path,
        path_kpoints=kpoint_path,
        path_energies=energy_path,
        path_metadata=path_metadata,
        metadata=metadata_path,
    )
end

end
