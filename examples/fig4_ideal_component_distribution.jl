using Serialization

using TMBGMagneticHybridWannier
using TMBGMagneticHybridWannier.TMBGFigure4Projection

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    cache_path = length(ARGS) >= 1 ? abspath(ARGS[1]) :
        joinpath(root, "results", "paper_q20_converged", "magnetic_hw.jls")
    output_directory = length(ARGS) >= 2 ? abspath(ARGS[2]) :
        joinpath(root, "results", "paper_q20_converged", "fig4")
    isfile(cache_path) || error(
        "Magnetic cache not found. Run examples/paper_example.jl first, " *
        "or pass its magnetic_hw.jls path as the first argument.",
    )
    result = deserialize(cache_path)
    outputs = write_figure4_projection(result, output_directory)
    println("grid=$(outputs.grid)")
    println("sources=$(outputs.sources)")
    println("metadata=$(outputs.metadata)")
    println(read(outputs.metadata, String))
end

main()
