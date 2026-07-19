using Serialization

using TMBGMagneticHybridWannier
using TMBGMagneticHybridWannier.TMBGMagneticHW
using TMBGMagneticHybridWannier.TMBGIdealComponent
using TMBGMagneticHybridWannier.TMBGIdealComponentProjection
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.TMBGZeroField
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.WangContinuum

const CASE_NAME = "tmbg_pq1_20_lg11_smax5"
const RECIPROCAL_GRID_SIZE = 11
const MAGNETIC_TRANSLATION_CUTOFF = 5
const ACTIVE_BAND_COUNT = 6

function build_tmbg_magnetic_model()
    parameters = TMBGParams(
        1.04,            # twist angle (degrees)
        3,               # source continuum-model shell cutoff
        0.70,            # kappa = w0/w1
        1,               # valley
        -30.0,           # hBN sublattice potential (meV)
        (1, 0, 0),       # hBN acts on the monolayer
        50.0,            # displacement potential (meV)
    )
    continuum_model = build_tmbg_wang(
        parameters;
        lg=RECIPROCAL_GRID_SIZE,
        periodic_G=false,
    )
    return build_magnetic_hw(
        continuum_model,
        FluxSpec(1, 20, 1);
        mesh_factor=1,
        nactive=ACTIVE_BAND_COUNT,
        nvec=collect(0:1),
        svec=collect(
            -MAGNETIC_TRANSLATION_CUTOFF:MAGNETIC_TRANSLATION_CUTOFF,
        ),
        gc=4,
    )
end

function validate_cached_model(result::MagneticHWModel)
    result.flux.p == 1 || error("Cached flux numerator is not p=1")
    result.flux.q == 20 || error("Cached flux denominator is not q=20")
    result.model.lg == RECIPROCAL_GRID_SIZE ||
        error("Cached reciprocal grid is not lg=$RECIPROCAL_GRID_SIZE")
    size(result.hw.WLS, 2) == ACTIVE_BAND_COUNT ||
        error("Cached active-band count is not $ACTIVE_BAND_COUNT")
    result.svec == collect(
        -MAGNETIC_TRANSLATION_CUTOFF:MAGNETIC_TRANSLATION_CUTOFF,
    ) || error(
        "Cached magnetic-translation cutoff is not " *
        "smax=$MAGNETIC_TRANSLATION_CUTOFF",
    )
    return result
end

function write_case_parameters(path::String)
    open(path, "w") do io
        println(io, "system=tMBG")
        println(io, "theta_deg=1.04")
        println(io, "kappa=0.70")
        println(io, "displacement_potential_meV=50.0")
        println(io, "hBN_sublattice_potential_meV=-30.0")
        println(io, "flux_p=1")
        println(io, "flux_q=20")
        println(io, "reciprocal_grid_lg=$RECIPROCAL_GRID_SIZE")
        println(io, "magnetic_translation_smax=$MAGNETIC_TRANSLATION_CUTOFF")
        println(io, "active_band_count=$ACTIVE_BAND_COUNT")
    end
    return path
end

function write_magnetic_spectrum(result::MagneticHWModel, path::String)
    open(path, "w") do io
        println(io, "k1,k2,magnetic_band,energy_meV")
        for i2 in 1:result.l2, i1 in 1:result.l1
            momentum_index = i1 + (i2 - 1) * result.l1
            for band in axes(result.spectrum, 1)
                println(io, join((
                    (i1 - 1) / result.l1,
                    (i2 - 1) / (result.l2 * result.flux.q),
                    band,
                    result.spectrum[band, momentum_index],
                ), ','))
            end
        end
    end
    return path
end

function main()
    repository_root = normpath(joinpath(@__DIR__, ".."))
    output_directory = isempty(ARGS) ?
        joinpath(repository_root, "results", CASE_NAME) : abspath(ARGS[1])
    mkpath(output_directory)

    parameter_path = write_case_parameters(
        joinpath(output_directory, "case_parameters.txt"),
    )
    cache_path = joinpath(output_directory, "magnetic_hw_cache.jls")
    result = if isfile(cache_path)
        println("Reusing magnetic hybrid-Wannier cache: $cache_path")
        validate_cached_model(deserialize(cache_path))
    else
        println("Building the converged tMBG magnetic hybrid-Wannier model")
        built = build_tmbg_magnetic_model()
        serialize(cache_path, built)
        built
    end

    spectrum_path = write_magnetic_spectrum(
        result,
        joinpath(output_directory, "magnetic_spectrum.csv"),
    )
    target_band = select_isolated_target_band(result)
    geometry = compute_ideal_component(
        result;
        target_band=target_band,
        correction_bands=[target_band - 1],
        reciprocal_shift_periodic=false,
    )
    trace_summary_path, geometry_path = write_trace_condition_results(
        geometry,
        result,
        output_directory,
    )
    overlap_summary_path, overlap_path = write_overlap_diagnostics(
        geometry,
        result,
        output_directory,
    )
    projection = write_ideal_component_projection(
        result,
        joinpath(output_directory, "ideal_component_projection");
        target_band=target_band,
    )

    println("parameters=$parameter_path")
    println("magnetic_spectrum=$spectrum_path")
    println("trace_condition_summary=$trace_summary_path")
    println("quantum_geometry=$geometry_path")
    println("overlap_summary=$overlap_summary_path")
    println("overlap_diagnostics=$overlap_path")
    println("momentum_resolved_projection=$(projection.momentum_weights)")
    println("projection_summary=$(projection.summary)")
    println(read(trace_summary_path, String))
end

main()
