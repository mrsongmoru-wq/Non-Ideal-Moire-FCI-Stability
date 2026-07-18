using Serialization

using TMBGMagneticHybridWannier
using TMBGMagneticHybridWannier.TMBGMagneticHW
using TMBGMagneticHybridWannier.TMBGIdealComponent
using TMBGMagneticHybridWannier.TMBGFigure4Projection
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.TMBGZeroField
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.WangContinuum

function option(name::String, default::String)
    prefix = "--$name="
    match = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(match) ? default : ARGS[match][(length(prefix) + 1):end]
end

integer_option(name::String, default::Int) = parse(Int, option(name, string(default)))

function write_spectrum(result::MagneticHWModel, path::String)
    open(path, "w") do io
        println(io, "k1,k2,band,energy_meV")
        for i2 in 1:result.l2, i1 in 1:result.l1
            ik = i1 + (i2 - 1) * result.l1
            for band in axes(result.spectrum, 1)
                println(io, join((
                    (i1 - 1) / result.l1,
                    (i2 - 1) / (result.l2 * result.flux.q),
                    band,
                    result.spectrum[band, ik],
                ), ','))
            end
        end
    end
end

function build_paper_model(
    ;
    quick::Bool=false,
    lg::Int=9,
    smax::Int=5,
    nactive::Int=6,
)
    params = TMBGParams(
        1.04,            # twist angle (degrees)
        3,               # hexagonal cutoff of the source continuum model
        0.70,            # kappa = w0/w1
        1,               # valley
        -30.0,           # hBN sublattice potential (meV)
        (1, 0, 0),       # hBN acts on the monolayer
        50.0,            # displacement potential (meV)
    )
    lg = quick ? 5 : lg
    q = quick ? 6 : 20
    smax = quick ? 2 : smax
    model = build_tmbg_wang(params; lg=lg, periodic_G=false)
    return build_magnetic_hw(
        model,
        FluxSpec(1, q, 1);
        mesh_factor=1,
        nactive=nactive,
        nvec=collect(0:1),
        svec=collect(-smax:smax),
        gc=4,
    )
end

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    quick = "--quick" in ARGS
    lg = integer_option("lg", 9)
    smax = integer_option("smax", 5)
    nactive = integer_option("nactive", 6)
    default_output = joinpath(
        root,
        "results",
        quick ? "quick" : "paper_q20_converged",
    )
    output_directory = abspath(option("output", default_output))
    mkpath(output_directory)
    cache_path = abspath(option("cache", joinpath(output_directory, "magnetic_hw.jls")))

    result = if isfile(cache_path)
        println("Reusing magnetic hWF cache: $cache_path")
        deserialize(cache_path)
    else
        println("Building magnetic hWF spectrum")
        built = build_paper_model(
            ;
            quick=quick,
            lg=lg,
            smax=smax,
            nactive=nactive,
        )
        serialize(cache_path, built)
        built
    end

    spectrum_path = joinpath(output_directory, "magnetic_spectrum.csv")
    write_spectrum(result, spectrum_path)

    target = select_paper_target_band(result)
    geometry = compute_ideal_component(
        result;
        target_band=target,
        correction_bands=[target - 1],
        reciprocal_shift_periodic=false,
    )
    summary_path, geometry_path =
        write_ideal_component(geometry, result, output_directory)
    overlap_summary_path, overlap_path =
        write_overlap_distribution(geometry, result, output_directory)
    figure4 = write_figure4_projection(
        result,
        joinpath(output_directory, "fig4");
        target_band=target,
    )

    println("spectrum=$spectrum_path")
    println("summary=$summary_path")
    println("geometry=$geometry_path")
    println("overlap_summary=$overlap_summary_path")
    println("overlaps=$overlap_path")
    println("fig4_grid=$(figure4.grid)")
    println("fig4_metadata=$(figure4.metadata)")
    println(read(summary_path, String))
end

main()
